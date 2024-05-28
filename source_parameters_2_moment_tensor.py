import numpy as np


def NED2ENU(aMomentTensor):
    """convert RTP moment tensor to NED
    derived from source code of obspy
    https://github.com/obspy/obspy/blob/master/obspy/imaging/source.py#L115
    """
    signs = [1, 1, 1, 1, -1, -1]
    indices = [1, 0, 2, 3, 5, 4]
    return np.array([sign * aMomentTensor[ind] for sign, ind in zip(signs, indices)]).T


def my_print(x):
    """print numpy array without []"""
    for row in x:
        print(" ".join(map(lambda x: "{:.8e}\t".format(x), row)))


def print_MT_src_50(MT):
    """print moment tensor in the form needed by SeiSsol source type 50
    MT: 6 components of moment tensor
    """
    a = -np.array([MT[k] for k in [0, 3, 4, 3, 1, 5, 4, 5, 2]]).reshape((3, 3))
    my_print(a)


class source:
    def __init__(self, strike, dip, rake, Garea, slip):
        self.strike = np.radians(strike)
        self.dip = np.radians(dip)
        self.rake = np.radians(rake)
        self.Garea = Garea
        self.slip = slip

    def compute_moment_tensor_NED(self):
        """compute equivalent moment tensor of given source parameter
        in NED convention (North East Down)
        """
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

        MomentTensor = np.zeros((6))
        # 0   1  2  3  4  5
        # xx,yy,zz,xy,xz,yz
        # http://gfzpublic.gfz-potsdam.de/pubman/item/escidoc:65580/component/escidoc:65579/IS_3.8_rev1.pdf (eq 5)
        # with x y z : NED
        MomentTensor[0] = -M0 * (sd * cl * s2s + s2d * sl * np.power(ss, 2))
        MomentTensor[1] = M0 * (sd * cl * s2s - s2d * sl * np.power(cs, 2))
        MomentTensor[2] = M0 * (s2d * sl)
        MomentTensor[3] = M0 * (sd * cl * c2s + 0.5 * s2d * sl * s2s)
        MomentTensor[4] = -M0 * (cd * cl * cs + c2d * sl * ss)
        MomentTensor[5] = -M0 * (cd * cl * ss - c2d * sl * cs)
        self.MomentTensor_NED = MomentTensor
        print("NN, EE, DD, NE, ND, ED")
        print(self.MomentTensor_NED)
        print_MT_src_50(self.MomentTensor_NED)

    def compute_moment_tensor_ENU(self):
        self.compute_moment_tensor_NED()
        self.MomentTensor_ENU = NED2ENU(self.MomentTensor_NED)
        print("EE, NN, UU, EN, EU, NU")
        print(self.MomentTensor_ENU)
        print_MT_src_50(self.MomentTensor_ENU)


src = source(strike=20, dip=70, rake=30, Garea=1e18, slip=1)
src.compute_moment_tensor_ENU()
