"""
Postprocessing of results
"""
import pickle
import os
import matplotlib.pyplot as plt
import numpy as np

from utils.utils import export_to, export_figure


class MSAPostprocessor:
    def __init__(self, path, flag3d=True, export=True):
        """
        :param path: str                    IDA results directory
        :param flag3d: bool                 3D modelling or not
        :param export: bool                 Exporting postprocessed IDA results or not
        """
        self.path = path
        self.flag3d = flag3d
        self.export = export

        # Alternative font size
        self.FONTSIZE = 10
        # Default color patterns
        self.color_grid = ['#840d81', '#6c4ba6', '#407bc1', '#18b5d8', '#01e9f5',
                           '#cef19d', '#a6dba7', '#77bd98', '#398684', '#094869']
        self.grayscale = ['#111111', '#222222', '#333333', '#444444', '#555555',
                          '#656565', '#767676', '#878787', '#989898', '#a9a9a9']

    def get_ground_motion_batches(self):
        return next(os.walk(self.path))[1]

    @staticmethod
    def read_pickle(path):
        with open(path, 'rb') as file:
            data = pickle.load(file)
        return data

    def msa(self, nst):
        if self.flag3d:
            n_directions = 2
        else:
            n_directions = 1

        gm_levels = self.get_ground_motion_batches()

        # outputs
        out = {}

        # For each level of excitation
        for level in gm_levels:
            try:
                mafe = float(level.split("-")[1])
                rp = int(1 / mafe)
            except:
                rp = int(level)

            out[rp] = {"1": {"acc": {}, "disp": {}, "drift": {}},
                       "2": {"acc": {}, "disp": {}, "drift": {}}}

            print(f"[LEVEL] {level}, Return period {rp} years")

            # For each direction
            for d in range(n_directions):
                for st in range(nst+1):
                    out[rp][str(d + 1)]["acc"][st] = []
                    out[rp][str(d + 1)]["disp"][st] = []
                    if st != nst:
                        out[rp][str(d + 1)]["drift"][st] = []

            # for each record level
            for record in next(os.walk(self.path / level))[-1]:
                if record.endswith("_part.pickle"):
                    continue
                """
                0 - accelerations [g], 1 - displacements [m], 2 - drifts
                each has a shape of d x s x r, where 
                    d stands for number of directions
                    s stands for number of storeys and floors
                    r stands for number of steps of the record
                """
                data = self.read_pickle(self.path / level / record)

                for d in range(n_directions):
                    for st in range(nst+1):
                        out[rp][str(d + 1)]["acc"][st].append(abs(max(data[0][d][st], key=abs)))
                        out[rp][str(d + 1)]["disp"][st].append(abs(max(data[1][d][st], key=abs)))
                        if st != nst:
                            out[rp][str(d + 1)]["drift"][st].append(abs(max(data[2][d][st], key=abs)))

        # Exporting
        if self.export:
            export_to(self.path / "outputs", out, "json")

        return out

    @staticmethod
    def get_return_periods(out):
        rp = []
        for period in out.keys():
            rp.append(int(period))
        rp.sort()

        return rp

    @staticmethod
    def get_edp(out, direction, storey, edptype, rps, factor=1.0):
        edp = []
        # For each return period
        for rp in rps:
            rp = str(rp)
            if edptype == "drift":
                edp.append(out[rp][str(direction+1)][edptype][str(storey)])
            else:
                # acc
                try:
                    critical = np.array(max(out[rp][str(1)][edptype][str(storey)], out[rp][str(2)][edptype][str(storey)])) \
                               * factor
                except:
                    critical = np.array(out[rp][str(1)][edptype][str(storey)]) * factor

                edp.append(list(critical))

        return edp

    @staticmethod
    def compute_pga(out, direction, storey, edptype, rps):
        pass

    def plot_vs_rp(self, edp, rp, filename, xlabel=None, ylabel=None, xlimit=None, ylimit=None, pflag=True):
        def median(lst):
            n = len(lst)
            s = sorted(lst)
            return (s[n // 2 - 1] / 2.0 + s[n // 2] / 2.0, s[n // 2])[n % 2] if n else None

        average = []
        for i in edp:
            average.append(median(i))

        fig, ax = plt.subplots(figsize=(4, 3), dpi=100)
        for i in range(len(edp)):
            x = edp[i]
            y = [rp[i]] * len(edp[i])
            if i == 0:
                label = "Records"
            else:
                label = None
            plt.scatter(x, y, color=self.color_grid[3], marker="o", label=label)

        plt.plot(average, rp, color=self.color_grid[0], marker="o", label="Median")
        plt.xlabel(xlabel, fontsize=self.FONTSIZE)
        plt.ylabel(ylabel, fontsize=self.FONTSIZE)
        plt.grid(True, which="major", axis='both', ls="--", lw=1.0)
        plt.grid(True, which="minor", axis='both', ls="--", lw=0.5)
        plt.yscale("log")
        if xlimit is not None and ylimit is not None:
            plt.xlim([0, xlimit])
            plt.ylim([1, ylimit])
        plt.rc('xtick', labelsize=self.FONTSIZE)
        plt.rc('ytick', labelsize=self.FONTSIZE)

        ax.legend(frameon=False, loc='best', fontsize=8)

        if pflag:
            plt.show()
        if filename:
            export_figure(fig, filename=filename, filetype="svg")

        plt.close()
