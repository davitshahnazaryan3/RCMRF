import matplotlib
import matplotlib.pyplot as plt
import pickle
from pathlib import Path

from utils_plotter import *


class Visualize:
    def __init__(self, export=False, filetype="emf", export_dir=None, flag=True):
        # Default color patterns
        self.color_grid = ['#840d81', '#6c4ba6', '#407bc1', '#18b5d8', '#01e9f5',
                           '#cef19d', '#a6dba7', '#77bd98', '#398684', '#094869']
        self.grayscale = ['#111111', '#222222', '#333333', '#444444', '#555555',
                          '#656565', '#767676', '#878787', '#989898', '#a9a9a9']

        # Set default font style
        font = {'size': 10}
        matplotlib.rc('font', **font)

        # Alternative font size
        self.FONTSIZE = 10

        # Exporting figures
        self.export = export
        self.filetype = filetype
        self.export_dir = export_dir
        self.flag = flag

    def plot_spo(self, filename, site="", direction=""):
        """
        SPO plotter
        :param filename: str
        :return: None
        """
        with open(filename, "rb") as f:
            spo = pickle.load(f)

        fig, ax = plt.subplots(figsize=(4, 3), dpi=100)
        plt.plot(spo[0] * 100, spo[1], color=self.grayscale[0])

        plt.xlabel("Top displacement [cm] ", fontsize=self.FONTSIZE)
        plt.ylabel('Base shear [kN]', fontsize=self.FONTSIZE)
        plt.grid(True, which="major", axis='both', ls="--", lw=1.0)
        plt.grid(True, which="minor", axis='both', ls="--", lw=0.5)
        plt.xlim([0, int(max(spo[0]) * 100) + 20])
        plt.ylim([0, int(max(spo[1])) + 300])
        plt.rc('xtick', labelsize=self.FONTSIZE)
        plt.rc('ytick', labelsize=self.FONTSIZE)
        # plt.legend(frameon=False, loc='upper right', fontsize=self.FONTSIZE)

        if self.flag:
            plt.show()

        if self.export:
            export_figure(fig, filename=self.export_dir / f"spo_{site}_{direction}", filetype=self.filetype)


if __name__ == "__main__":
    path = Path.cwd().parents[0]
    export_dir = path / "sample/RCMRF/figs"
    create_folder(export_dir)

    site = "low"
    direction = "2"
    spo_model = path / f"sample/RCMRF/SPO_{site}_{direction}.pickle"

    viz = Visualize(export=True, filetype="png", export_dir=export_dir, flag=True)
    viz.plot_spo(spo_model, site=site, direction=direction)
