import json
import pickle

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from postprocess.utils_plotter import *


class Visualize:
    def __init__(self, export=False, filetype="emf", export_dir=None, flag=True):
        """
        Initialize
        :param export: bool             Exporting to file?
        :param filetype: str            File extension
        :param export_dir: Path         Path to export to
        :param flag: bool               Show figure?
        """
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

    @staticmethod
    def read_file(filename):

        with open(filename, "rb") as f:
            if filename.suffix == ".json":
                # json format
                data = json.load(f)
            else:
                # pickle format
                data = pickle.load(f)

        return data

    def plot_spo(self, filename, name=None, labels=None):
        """
        SPO plotter
        :param filename: Path
        :param name: str
        :param labels: List[str]
        :return: None
        """
        # Initialize figure
        fig, ax = plt.subplots(figsize=(4, 3), dpi=100)

        if not isinstance(filename, list):
            # A single graph
            spo = self.read_file(filename)

            maxvalx = max(spo[0]) * 100
            maxvaly = max(spo[1])
            plt.plot(np.array(spo[0]) * 100, spo[1], color=self.grayscale[0])

        else:
            # multiple graphs
            cnt = 0
            maxvalx = -float("inf")
            maxvaly = -float("inf")
            for i in range(len(filename)):
                model = filename[i]
                if labels is not None:
                    label = labels[i]
                else:
                    label = None

                spo = self.read_file(model)

                plt.plot(np.array(spo[0]) * 100, spo[1], color=self.color_grid[cnt], label=label)
                cnt += 2

                if maxvalx < max(spo[0]) * 100:
                    maxvalx = max(spo[0]) * 100
                if maxvaly < max(spo[1]):
                    maxvaly = max(spo[1])

        plt.xlabel("Top displacement [cm]", fontsize=self.FONTSIZE)
        plt.ylabel('Base shear [kN]', fontsize=self.FONTSIZE)
        plt.grid(True, which="major", axis='both', ls="--", lw=1.0)
        plt.grid(True, which="minor", axis='both', ls="--", lw=0.5)
        plt.xlim([0, int(maxvalx) + 20])
        plt.ylim([0, int(maxvaly) + 300])
        plt.rc('xtick', labelsize=self.FONTSIZE)
        plt.rc('ytick', labelsize=self.FONTSIZE)
        plt.legend(frameon=False, loc='upper right', fontsize=self.FONTSIZE)

        if self.flag:
            plt.show()

        if self.export:
            try:
                # If inkscape exists
                export_figure(fig, filename=self.export_dir / f"spo_{name}", filetype=self.filetype)

            except:
                export_figure_basic(fig, filename=self.export_dir / f"spo_{name}", filetype=self.filetype)
