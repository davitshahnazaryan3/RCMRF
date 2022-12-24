"""
user defines input arguments
Main file to run the software
"""
import numpy as np


class Input:
    # Initialize input arguments
    heights = [4.5, 3.5, 3.5, 3.5]
    nst = len(heights)
    spans_x = [6., 4.5, 3.5, 6., 3.5, 4.5, 6.]
    n_bays = len(spans_x)
    spans_y = [6., 4., 6.]
    fy = 435.
    elastic_modulus_steel = 200000.
    eps_y = fy / elastic_modulus_steel
    fc = 25.
    Ec = (3320 * fc**0.5 + 6900) * 1000.
    n_seismic = None
    n_gravity = None
    configuration = "space"

    site = None
    site_idx = 0

    # Live load
    live_load = 2.

    # Dead load
    dead_load = .2 * 25. + 1.

    # Snow loads, Milano/AnconaSaAvg -> 1.5, L'Aquila -> 2.2
    snow_load = [1.5, 1.5, 2.2]

    inputs = {"loads": [], "seismic": []}

    masses = None

    def __init__(self, site, flag3d=True):
        """
        initializes the input functions
        """
        self.flag3d = flag3d
        self.site = site

        if site.lower() == "milano":
            self.site_idx = 0
        elif site.lower() == "ancona":
            self.site_idx = 1
        elif site.lower() == "laquila":
            self.site_idx = 2
        else:
            raise ValueError("Wrong site location!")

        # Important only for 2D modelling
        if self.configuration == "perimeter" or not flag3d:
            # Masses will be subdivided between two seismic frames
            self.n_seismic = 2
            self.n_gravity = int(len(self.spans_y) - 1)
        else:
            # Masses will be considered for the entirety of the building considering all seismic frames
            self.n_gravity = 0
            self.n_seismic = 1

    def get_masses(self):
        # Floor area
        area = sum(self.spans_x) * sum(self.spans_y)
        self.masses = np.zeros(4)
        for st in range(self.nst):
            self.masses[st] = area / 9.81 * self.inputs["seismic"][st]

        return self.masses
