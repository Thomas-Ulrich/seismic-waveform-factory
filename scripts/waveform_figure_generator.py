import os

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.offsetbox import AnchoredText
from obspy import Trace
from obspy.signal.cross_correlation import correlate, xcorr_max
from obspy.signal.tf_misfit import eg, pg

plt.rcParams["font.family"] = "sans"


def remove_top_right_axes(ax):
    # remove top right axis
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()


def autoscale_y(ax, margin=0.1):
    """This function rescales the y-axis based on the data that is visible given the
    current xlim of the axis.

    ax -- a matplotlib axes object
    margin -- the fraction of the total height of the y-data to pad the upper and
    lower ylims
    """

    def get_bottom_top(line):
        xd = line.get_xdata()
        yd = line.get_ydata()
        lo, hi = ax.get_xlim()
        y_displayed = yd[((xd > lo) & (xd < hi))]
        if len(y_displayed) > 0:
            h = np.max(y_displayed) - np.min(y_displayed)
            bot = np.min(y_displayed) - margin * h
            top = np.max(y_displayed) + margin * h
        else:
            bot, top = np.inf, -np.inf
        return bot, top

    lines = ax.get_lines()
    bot, top = np.inf, -np.inf

    for line in lines:
        new_bot, new_top = get_bottom_top(line)
        if new_bot < bot:
            bot = new_bot
        if new_top > top:
            top = new_top
    if np.isfinite(bot) and np.isfinite(top):
        ax.set_ylim(bot, top)


def add_fault_components(streams, fault_strike):
    """
    Given a list of ObsPy streams and a fault strike angle,
    rotate N/E components into fault-parallel (f) and fault-normal (o) components
    and add them to each stream.

    Parameters
    ----------
    streams : list
        List of ObsPy Stream objects to process.
        (You can include st_obs here too.)
    fault_strike : float
        Fault strike angle in degrees.
    """
    theta_fp = np.deg2rad(fault_strike)
    theta_fn = np.deg2rad(fault_strike + 90)

    for st in streams:
        trN = st.select(component="N")[0]
        trE = st.select(component="E")[0]

        # Rotate
        fp = trN.data * np.cos(theta_fp) + trE.data * np.sin(theta_fp)
        fn = trN.data * np.cos(theta_fn) + trE.data * np.sin(theta_fn)

        # Create new traces
        trFP = Trace(data=fp, header=trN.stats)
        trFP.stats.channel = trN.stats.channel[:-1] + "f"  # fault parallel

        trFN = Trace(data=fn, header=trN.stats)
        trFN.stats.channel = trN.stats.channel[:-1] + "o"  # fault normal

        # Append to the same stream
        st += trFP
        st += trFN

    return streams


class WaveformFigureGenerator:
    def __init__(
        self,
        general_cfg,
        plt_cfg,
        n_kinematic_models,
    ):
        mapping = {"normal": "o", "parallel": "f"}
        self.components = [mapping.get(c, c) for c in plt_cfg["components"]]
        self.fault_strike = plt_cfg["fault_strike"]
        if "o" in self.components or "f" in self.components:
            assert self.fault_strike is not None, (
                "fault_strike parameter required for components='f' or 'o'"
                "(fault-parallel and normal)"
            )
        self.gen_cfg = general_cfg
        self.plt_cfg = plt_cfg
        self.plot_type = plt_cfg["type"]
        self.signal_kind = general_cfg["kind"]
        self.t_before = -plt_cfg["t_before"]
        self.t_after = plt_cfg["t_after"]
        self.filter_fmin = 1.0 / plt_cfg["filter_tmax"]
        self.filter_fmax = 1.0 / plt_cfg["filter_tmin"]
        self.enabled = plt_cfg["enabled"]
        self.ncol_per_component = plt_cfg["ncol_per_component"]
        self.ncomp = len(self.components)
        self.shift_match_correlation = plt_cfg["shift_match_correlation"]
        nstations = len(general_cfg["stations"])
        self.init_several_stations_figure(nstations)
        self.n_kinematic_models = n_kinematic_models
        self.kind_misfit = general_cfg["misfit"]
        self.init_gof_pandas_df()
        self.scaling = plt_cfg["scaling"]
        self.colors = general_cfg["line_colors"]
        self.line_widths = general_cfg["line_widths"]
        self.relative_offset = general_cfg["relative_offset"]
        self.annotations = plt_cfg["annotations"]
        self.global_legend_labels = general_cfg["global_legend_labels"]
        self.estimated_travel_time = 0.0

    def init_gof_pandas_df(self):
        columns = ["station", "distance", "azimuth"]
        if len(self.components) > 1:
            prefix = "".join(self.components)
            for i in range(self.n_kinematic_models):
                columns += [f"{self.plot_type}_{prefix}{i}"]
        for i in range(self.n_kinematic_models):
            columns += [f"{self.plot_type}_{comp}{i}" for comp in self.components]
        self.gof_df = pd.DataFrame(columns=columns)

    def set_estimated_travel_time(self, travel_time):
        self.estimated_travel_time = travel_time

    def init_several_stations_figure(self, nstations):
        nrow = int(np.ceil(nstations / self.ncol_per_component))
        ncol = self.ncol_per_component * self.ncomp
        # for comparing signals on the same figure
        if self.plot_type in ["p", "sh"]:
            height_one_plot = 1.7
            generic_plot = False
        else:
            height_one_plot = 1.5
            generic_plot = True
        fig, axarr = plt.subplots(
            nrow,
            ncol,
            figsize=(ncol * 4, nrow * height_one_plot),
            dpi=160,
            sharex=False,
            sharey=False,
            squeeze=False,
        )
        for j in range(ncol):
            for i in range(nrow):
                axi = axarr[i, j]
                if self.plt_cfg["normalize"]:
                    axi.set_yticks([])
                if j > 0 and generic_plot:
                    axi.set_yticks([])
                    axi.spines["left"].set_visible(False)
                    axi.get_shared_y_axes().joined(axi, axarr[i, 0])
                remove_top_right_axes(axi)
                axi.tick_params(axis="x", zorder=3)
                if i < nrow - 1 and generic_plot:
                    axi.spines["bottom"].set_visible(False)
                    axi.set_xticks([])
                if j * nrow + i >= self.ncomp * nstations:
                    axi.spines["left"].set_visible(False)
                    axi.spines["bottom"].set_visible(False)
                    axi.set_xticks([])
                    axi.set_yticks([])
        self.fig, self.axarr = fig, axarr

    def add_global_legend(self):
        # Add an invisible axis for the legend above the figure
        self.fig.subplots_adjust(top=0.88)  # Make some room at the top
        self.legend_ax = self.fig.add_axes([0.1, 0.9, 0.8, 0.05], frameon=False)
        self.legend_ax.axis("off")
        # Define one Line2D for each synthetic model and observation
        handles = []
        nlabels = len(self.global_legend_labels)
        ncol = self.ncol_per_component * self.ncomp
        for idx in range(nlabels):
            if idx > len(self.line_widths) - 1:
                color = "k"
                line_width = self.line_widths[-1]
            else:
                color = color = self.colors[idx]
                line_width = self.line_widths[idx]
            line = mlines.Line2D([], [], color=color, linewidth=line_width)
            handles.append(line)

        self.legend_ax.legend(
            handles,
            self.global_legend_labels,
            loc="center",
            ncol=ncol,
            frameon=False,
        )

    def compute_max_abs_value_trace(self, trace, reftime):
        max_abs_value = float("-inf")
        times = trace.times(reftime=reftime)
        mask = (times >= self.t_before) & (times <= self.t_after)
        filtered_data = trace.data[mask]
        if filtered_data.size > 0:
            max_abs_value = max(np.max(filtered_data), -np.min(filtered_data))
        return max_abs_value

    def compute_max_abs_value(self, streams, reftime):
        max_abs_value = float("-inf")
        for stream in streams:
            for comp in self.components:
                traces = stream.select(component=comp)
                if not traces:
                    continue
                trace = traces[0]
                max_trace = self.compute_max_abs_value_trace(trace, reftime)
                max_abs_value = max(max_abs_value, max_trace)
        if max_abs_value == float("-inf"):
            print(
                "No data points found in the specified time range for the"
                " specified components."
            )
            return 0.0
        else:
            return max_abs_value

    def compute_scaling(self, trace, reftime):
        annot = ""
        if self.plt_cfg["normalize"]:
            stmax = self.compute_max_abs_value_trace(trace, reftime)
            annot = f"{stmax:.2}"
            if stmax <= 0:
                stmax = 1.0
            scaling = 1 / stmax if stmax != 0.0 else 1.0
        elif self.plot_type in ["p", "sh"]:
            # micro meters (/s if velocity)
            scaling = 1e6
        else:
            scaling = 1
        scaling *= self.scaling
        return scaling, annot

    def add_plot_station(self, st_obs0, lst, reftime, ista):
        network = st_obs0[0].stats.network
        station = st_obs0[0].stats.station
        st_obs = st_obs0.copy()
        lst_copy = [st.copy() for st in lst]

        if "T" in self.components:
            for st in [*lst_copy, st_obs]:
                st.rotate(method="NE->RT")

        if "f" in self.components or "o" in self.components:
            streams = [*lst_copy, st_obs]
            add_fault_components(streams, self.fault_strike)

        st_obs = st_obs.split()
        for myst in [*lst_copy, st_obs]:
            myst.detrend("demean")
            myst.detrend("linear")
            if self.plt_cfg["taper"]:
                myst.taper(max_percentage=0.05, type="hann")
            myst.filter(
                "bandpass",
                freqmin=self.filter_fmin,
                freqmax=self.filter_fmax,
                corners=4,
                zerophase=True,
            )
        st_obs.merge()
        if self.plt_cfg["normalize"]:
            offset = self.relative_offset
        else:
            offset = (
                self.compute_max_abs_value([*lst_copy, st_obs], reftime)
                * self.relative_offset
            )

        nrows = self.axarr.shape[0]
        ins = ista % nrows

        def compute_j0(j):
            return j + ista // nrows

        ylabel = f"{network}.{station}"
        if self.plot_type == "generic":
            self.axarr[ins, 0].set_ylabel(ylabel)

        for j, comp in enumerate(self.components):
            j0 = compute_j0(j)
            if self.plot_type in ["p", "sh"]:
                self.axarr[ins, j0].set_ylabel(ylabel)
            vmax_annot = []
            nst = len(lst_copy) + 1
            # initialize ist if no synthetics
            ist = 0
            for ist, st in enumerate(lst_copy):
                strace = st.select(component=comp)[0]
                scaling, annot = self.compute_scaling(strace, reftime)
                shift_s = self.compute_shift_correlation(
                    st, st_obs, self.components, reftime
                )
                vmax_annot.append(annot)

                self.axarr[ins, j0].plot(
                    strace.times(reftime=reftime + shift_s),
                    scaling * strace.data + (nst - ist - 1) * offset,
                    self.colors[ist],
                    linewidth=self.line_widths[ist],
                )
            otraces = st_obs.select(component=comp)
            if not otraces:
                otrace = st_obs[0].copy()
                otrace.data *= np.nan
            else:
                otrace = otraces[0]
            scaling, annot = self.compute_scaling(otrace, reftime)
            vmax_annot.append(annot)
            self.axarr[ins, j0].plot(
                otrace.times(reftime=reftime),
                scaling * otrace.data,
                "k",
                linewidth=self.line_widths[ist],
            )
            self.axarr[ins, j0].set_xlim([self.t_before, self.t_after])
            # rescale to view range
            if len(self.components) == 1:
                autoscale_y(self.axarr[ins, j0])

        # Compute rMRS misfit and print it on plot
        if self.annotations["distance_unit"] == "degree":
            dist = otrace.stats.distance
            distance_unit = "°"
        else:
            dist = otrace.stats.distance_km
            distance_unit = "km"

        azimuth = otrace.stats.back_azimuth
        n_kinematic_models = len(lst_copy)
        temp_dic = {
            "station": station,
            "distance": dist,
            "azimuth": azimuth,
            "distance_unit": distance_unit,
        }
        for j, comp in enumerate(self.components):
            j0 = compute_j0(j)
            ymin0, ymax0 = self.axarr[ins, j0].get_ylim()
            gofstrings = []
            annot = []

            for ist, myst in enumerate(lst_copy):
                try:
                    gof, y0 = self.compute_misfit(myst, st_obs, comp, reftime)
                except ValueError:
                    gof, y0 = 0, 0
                temp_dic[f"{self.plot_type}_{comp}{ist}"] = gof
                if "misfit" in self.annotations["fields"]:
                    gofstrings += [f"{gof:.2f}"]
                if j == 0 and ist == n_kinematic_models - 1:
                    if "distance" in self.annotations["fields"]:
                        annot += [f"d:{dist:.0f}{distance_unit}"]
                    if "azimuth" in self.annotations["fields"]:
                        annot += [f"a:{azimuth:.0f}°"]
                if ist == n_kinematic_models - 1:
                    if "misfit" in self.annotations["fields"]:
                        gofstring = "\n" + " ".join(gofstrings)
                        annot += [f"{gofstring}"]
                    if self.plt_cfg["normalize"]:
                        annot += ["\n".join(vmax_annot)]
                annotations = " ".join(annot)
                loc = "upper" if (y0 - ymin0) / (ymax0 - ymin0) < 0.5 else "lower"
                if len(annot):
                    anchored_text = AnchoredText(
                        annotations, loc=loc + " left", frameon=False
                    )
                    self.axarr[ins, j0].add_artist(anchored_text)

        # compute average gof if several components
        if len(self.components) > 1:
            prefix = "".join(self.components)
            for i in range(self.n_kinematic_models):
                lgof = [
                    temp_dic[f"{self.plot_type}_{comp}{i}"] for comp in self.components
                ]
                av = sum(lgof) / self.ncomp
                temp_dic[f"{self.plot_type}_{prefix}{i}"] = av

        self.gof_df.loc[len(self.gof_df)] = temp_dic

    def _prepare_traces(self, st, st_obs, comp, reftime):
        """
        Common preprocessing: select, trim, and interpolate traces for one component.
        Returns (strace, otrace, f0, shiftmax).
        If preparation fails, returns (None, None, None, None).
        """
        straces = st.select(component=comp)
        otraces = st_obs.select(component=comp)

        if not straces or not otraces:
            return None, None, None, None

        strace = straces[0].copy()
        otrace = otraces[0].copy()

        # Find common time window
        start_osTrace = max(strace.stats.starttime, otrace.stats.starttime)
        end_osTrace = min(strace.stats.endtime, otrace.stats.endtime)
        start_time_interp = max(reftime + self.t_before, start_osTrace)
        end_time_interp = min(reftime + self.t_after, end_osTrace)

        # Interpolation frequency
        f0 = self.filter_fmax * 4.0
        npts_interp = int(np.floor((end_time_interp - start_time_interp) * f0)) + 1
        if npts_interp < 2:
            return None, None, None, None

        # Interpolate both traces
        strace.interpolate(
            sampling_rate=f0, starttime=start_time_interp, npts=npts_interp
        )
        otrace = otrace.split()
        otrace.interpolate(
            sampling_rate=f0, starttime=start_time_interp, npts=npts_interp
        )
        otrace = otrace.merge()[0]

        # Max allowed shift
        shift_sec_max = (
            100.0
            if self.plot_type == "generic"
            else max(2.5, 0.025 * self.estimated_travel_time)
        )
        shiftmax = int(shift_sec_max * f0)

        return strace, otrace, f0, shiftmax

    def compute_shift_correlation(self, st, st_obs, components, reftime):
        shift = 0
        for comp in self.components:
            shift += self.compute_shift_correlation_one_comp(st, st_obs, comp, reftime)
        return shift / len(self.components)

    def compute_shift_correlation_one_comp(self, st, st_obs, comp, reftime):
        if not self.shift_match_correlation:
            return 0

        strace, otrace, f0, shiftmax = self._prepare_traces(st, st_obs, comp, reftime)
        if strace is None:
            return 0

        cc = correlate(strace, otrace, shift=shiftmax)
        shift, gof = xcorr_max(cc, abs_max=False)
        return shift / f0

    def compute_misfit(self, st, st_obs, comp, reftime):
        strace, otrace, f0, shiftmax = self._prepare_traces(st, st_obs, comp, reftime)
        if strace is None:
            # Build dummy otrace if missing
            otrace = st_obs[0].copy()
            otrace.data *= 0
            return 0, otrace.data[0] * self.scaling

        def nanrms(x, axis=None):
            return np.sqrt(np.nanmean(x**2, axis=axis))

        if self.kind_misfit == "normalized_rms":
            gof = nanrms(strace.data - otrace.data) / nanrms(otrace.data)
            # transform misfit to goodness of fit
            gof = np.exp(-gof)
        elif self.kind_misfit == "min_shifted_normalized_rms":
            gof = np.inf
            # best_shift = 0
            for shift in range(-shiftmax, shiftmax + 1):
                if shift >= 0:
                    s_data = strace.data[shift:]
                    o_data = otrace.data[:-shift] if shift > 0 else otrace.data
                else:
                    shift_abs = abs(shift)
                    s_data = strace.data[:-shift_abs]
                    o_data = otrace.data[shift_abs:]
                current_gof = nanrms(s_data - o_data) / nanrms(o_data)
                if current_gof < gof:
                    gof = current_gof
                    # best_shift = shift
            # transform misfit to goodness of fit
            gof = np.exp(-gof)
        elif self.kind_misfit == "cross-correlation":
            cc = correlate(strace, otrace, shift=shiftmax)
            shift, gof = xcorr_max(cc, abs_max=False)
        elif self.kind_misfit == "time-frequency":
            gof_envolope = eg(
                strace.data,
                otrace.data,
                1 / f0,
                fmin=self.filter_fmin,
                fmax=self.filter_fmax,
                nf=100,
                w0=6,
                norm="global",
                a=1.0,
                k=1.0,
            )
            gof_phase = pg(
                strace.data,
                otrace.data,
                1 / f0,
                fmin=self.filter_fmin,
                fmax=self.filter_fmax,
                nf=100,
                w0=6,
                norm="global",
                a=1.0,
                k=1.0,
            )
            gof = (2.0 * gof_envolope + gof_phase) / 3.0
        else:
            raise NotImplementedError(self.kind_misfit)
        return gof, otrace.data[0] * self.scaling

    def finalize_and_save_fig(self, fname=None):
        if fname is None:
            setup_name = self.gen_cfg["setup_name"]
            ext = self.gen_cfg["figure_extension"]
            fname = f"plots/{setup_name}_{self.signal_kind}_{self.plot_type}.{ext}"

        if self.global_legend_labels:
            self.add_global_legend()
        if self.plot_type == "generic":
            direction = {"E": "EW", "N": "NS", "Z": "UD", "f": "FP", "o": "FN"}
            for j, comp in enumerate(self.components):
                self.axarr[0, j].set_title(direction[comp])
            self.axarr[-1, -1].set_xlabel("Time (s)")
        elif self.plot_type == "sh":
            self.axarr[0, 0].set_title("T")
            self.axarr[-1, -1].set_xlabel("time relative to SH arrival (s)")
        elif self.plot_type == "p":
            direction = {"E": "EW", "N": "NS", "Z": "UD"}
            for j, comp in enumerate(self.components):
                self.axarr[0, j].set_title(direction[comp])
            self.axarr[-1, -1].set_xlabel("time relative to P arrival (s)")
        for j in range(self.ncol_per_component):
            self.fig.align_ylabels(self.axarr[:, j])
        self.fig.subplots_adjust(wspace=0.3 if self.ncol_per_component > 1 else 0.1)
        self.fig.savefig(fname, bbox_inches="tight")
        print(f"done generating {fname}")
        full_path = os.path.abspath(fname)
        print(f"full path: {full_path}")
