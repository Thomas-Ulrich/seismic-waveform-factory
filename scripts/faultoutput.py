import numpy as np
import seissolxdmf
import cmt


class FaultOutput:
    def __init__(self, fname, idt=None):
        self.sx = seissolxdmf.seissolxdmf(fname)
        self.xyz = self.sx.ReadGeometry()
        self.connect = self.sx.ReadConnect()
        self.nElements = self.connect.shape[0]
        if idt:
            self.ndt = idt[0] + 1
        else:
            self.ndt = self.sx.ndt
        if "fault-tag" in self.sx.ReadAvailableDataFields():
            self.fault_tags = self.sx.Read1dData(
                "fault-tag", self.nElements, isInt=True
            ).T
        else:
            self.fault_tags = np.zeros((self.nElements,), dtype=int) + 3
        self.unique_fault_tags = np.unique(self.fault_tags)

    def read_final_slip(self):
        """Store slip variable."""
        self.slip = self.sx.ReadData("ASl", self.ndt - 1)

    def get_rupture_time(self):
        """Return RT variable."""
        return self.sx.ReadData("RT", self.ndt - 1)

    def compute_strike_dip(self, refVector):
        """Compute strike and dip angles of each face of fault output."""
        print("computing strike and dip angles...")
        un = np.cross(
            self.xyz[self.connect[:, 1], :] - self.xyz[self.connect[:, 0], :],
            self.xyz[self.connect[:, 2], :] - self.xyz[self.connect[:, 0], :],
        )
        norm = np.apply_along_axis(np.linalg.norm, 1, un)
        un = un / norm.reshape((self.nElements, 1))
        un = un.T
        my_sign = np.sign(np.dot(-un.T, refVector))
        if np.any(my_sign == 0):
            raise ValueError("Wrong refVector: is locally normal to fault normal")
        un[:, :] = un[:, :] * my_sign
        # compute strike and dip direction
        us = np.zeros(un.shape)
        us[0, :] = -un[1, :]
        us[1, :] = un[0, :]
        norm = np.apply_along_axis(np.linalg.norm, 0, us)
        us = us / norm
        ud = np.cross(un.T, us.T).T

        self.strike = np.arctan2(us[0, :], us[1, :])
        idsNumericalError = np.where(ud[2, :] > 1)
        ud[2, idsNumericalError] = 1
        self.dip = np.arcsin(ud[2, :])

    def compute_rake(self, invertSld):
        """Compute rake angle of fault slip."""
        Sls = self.sx.ReadData("Sls", self.ndt - 1)
        Sld = self.sx.ReadData("Sld", self.ndt - 1)
        if invertSld:
            Sld = -Sld
        # the minus sign on Sls is due to SeisSol convention
        # of positive strike for  right-lateral faulting
        self.rake = np.arctan2(Sld, -Sls)

    def compute_face_area(self):
        """Compute area of each triangle."""
        print("computing area...")
        cross0 = np.cross(
            self.xyz[self.connect[:, 1], :] - self.xyz[self.connect[:, 0], :],
            self.xyz[self.connect[:, 2], :] - self.xyz[self.connect[:, 0], :],
        )
        self.face_area = 0.5 * np.apply_along_axis(np.linalg.norm, 1, cross0)

    def evaluate_G(self, mu_description):
        if mu_description.endswith(".txt"):
            print(f"loading 1D velocity file {mu_description}")
            depthmu = np.loadtxt(mu_description)
            depthmu = depthmu[depthmu[:, 0].argsort()]
            self.G = np.interp(self.xyzc[:, 2], depthmu[:, 0], depthmu[:, 1])

        elif mu_description.endswith(".yaml") or mu_description.endswith(".yml"):
            import easi

            regions = np.ones((self.nElements,))
            print(
                "Warning: assuming region 1 for all cells when evaluating mu with easi"
            )
            out = easi.evaluate_model(self.xyzc, regions, ["mu"], mu_description)
            self.G = out["mu"]
        else:
            self.G = float(mu_description)

    def compute_Garea(self, G):
        """Compute rigidity times area of each triangle."""
        self.Garea = G * self.compute_face_area()

    def compute_barycenter_coords(self):
        """Compute coordinates of fault element barycenter."""
        print("computing barycenters...")
        self.xyzc = (
            self.xyz[self.connect[:, 0], :]
            + self.xyz[self.connect[:, 1], :]
            + self.xyz[self.connect[:, 2], :]
        ) / 3.0

    def compute_face_moment_rate_from_ASl(self):
        """Compute moment rate function of each face element of the fault output by
        derivating accumulated slip output."""
        print("computing moment rate function from accumulated slip...")
        try:
            self.dt = self.sx.ReadTimeStep()
        except NameError:
            # there is only one sample in fault output (extract of the last time step)
            self.dt = 1.0
        self.face_moment_rate = np.zeros((self.ndt, self.nElements))
        # if too many elements this may generate a memory overflow, therefore the if
        if self.ndt * self.nElements < 10e6:
            ASl = self.sx.ReadData("ASl")
            self.face_moment_rate[1:, :] = np.diff(ASl, axis=0) / self.dt
        else:
            print("using a slower but more memory efficient implementation")
            slip1 = self.sx.ReadData("ASl", 0)
            for i in range(1, self.ndt):
                print(i, end=" ")
                slip0 = np.copy(slip1)
                slip1 = self.sx.ReadData("ASl", i)
                self.face_moment_rate[i, :] = (slip1 - slip0) / self.dt
            print()
        self.face_moment_rate = self.face_moment_rate * self.Garea

    def compute_face_moment_rate_from_slip_rate(self, fo_SR):
        """Compute moment rate function of each face elemnt of the fault output using
        directly the slip rate variables.

        high sampling required!!!
        """
        print("computing moment rate function from SR time histories")
        SRs = fo_SR.sx.ReadData("SRd")
        SRd = fo_SR.sx.ReadData("SRd")
        same_connect = np.allclose(fo_SR.connect, self.connect)
        same_xyz = np.allclose(fo_SR.xyz, self.xyz)
        if not (same_xyz and same_connect):
            raise ValueError(
                "main output file and SR file no have have the same \
                 geometry or connect"
            )
        self.dt = fo_SR.sx.ReadTimeStep()
        self.ndt = fo_SR.sx.ndt
        self.face_moment_rate = self.Garea * np.sqrt(SRs**2 + SRd**2)

    def compute_face_moment_tensor_NED(self):
        """Compute equivalent moment tensor of each face of the fault output in NED
        convention (North East Down)"""
        cs = np.cos(self.strike)
        c2s = np.cos(2.0 * self.strike)
        cd = np.cos(self.dip)
        c2d = np.cos(2.0 * self.dip)
        cl = np.cos(self.rake)

        ss = np.sin(self.strike)
        s2s = np.sin(2.0 * self.strike)
        sd = np.sin(self.dip)
        s2d = np.sin(2.0 * self.dip)
        sl = np.sin(self.rake)

        M0 = self.Garea * self.slip

        moment_tensor = np.zeros((6, self.nElements))
        # 0   1  2  3  4  5
        # xx,yy,zz,xy,xz,yz
        # http://gfzpublic.gfz-potsdam.de/pubman/item/escidoc:65580/component/
        # escidoc:65579/IS_3.8_rev1.pdf (eq 5)
        # with x y z : NED
        moment_tensor[0, :] = -M0 * (sd * cl * s2s + s2d * sl * np.power(ss, 2))
        moment_tensor[1, :] = M0 * (sd * cl * s2s - s2d * sl * np.power(cs, 2))
        moment_tensor[2, :] = M0 * (s2d * sl)
        moment_tensor[3, :] = M0 * (sd * cl * c2s + 0.5 * s2d * sl * s2s)
        moment_tensor[4, :] = -M0 * (cd * cl * cs + c2d * sl * ss)
        moment_tensor[5, :] = -M0 * (cd * cl * ss - c2d * sl * cs)
        self.face_moment_tensor = moment_tensor

    def compute_equivalent_point_source_subfault(self, ids):
        """Compute properties of the equivalent moment tensor of a subset ids of faces
        of the fault outputs."""
        moment_rate = np.sum(self.face_moment_rate[:, ids], axis=1)
        # Note we do not use np.trapz here because we just want to revert our derivation
        mom = np.sum(moment_rate) * self.dt
        if abs(mom) < 1e-3:
            return None, None, None, None
        norm_moment_rate = moment_rate / mom

        moment_tensor = cmt.compute_moment_tensor(self.face_moment_tensor[:, ids])
        m0all = cmt.compute_seismic_moment(moment_tensor)

        # Correct the moment tensor to the get actual seismic moment
        # In fact, m0all can be smaller than mom:
        # For instance if 2 faces of similar orientation have opposite fault slip vector
        # mom will be 2*M0 and m0all will be zero
        if m0all > 0:
            corrected_moment_tensor = (mom / m0all) * moment_tensor

        face_moment = np.sum(self.face_moment_rate[:, ids], axis=0) * self.dt
        xyzc = np.average(self.xyzc[ids, :], axis=0, weights=face_moment)
        return mom, norm_moment_rate, corrected_moment_tensor, xyzc
