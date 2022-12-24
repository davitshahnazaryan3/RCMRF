"""
user defines input arguments
Main file to run the software
"""
import numpy as np


class Input:
    """
    * Those are the only inputs to be updated for your Model scenario:
    heights: list
        Heights of the buildings in m
    spans_x: list
        Spans in X direction in m
    spans_y: list
        Spans in Y direction in m (optional, required only for 3D analysis)

    -- Material properties
    fy: float
        Yield strength of reinforcement
    elastic_modulus_steel: float
        Elastic modulus of reinforcement
    fc: float
        Compressive strength of concrete
    -- Gravity loads (factored)
    live_load: float
    dead_load: float
    roof_live: float
    """
    # Initialize input arguments
    heights = [4.0, 3.0, 3.0, 3.0, 3.0, 3.0]
    spans_x = [3.5, 4.5, 3.5, 4., 3.5, 4.5, 3.5]
    fy = 435.
    elastic_modulus_steel = 200000.
    fc = 25.
    live_load = 2.
    dead_load = .15 * 25. + 1.
    roof_load = 2.2

    # Computed inputs
    spans_y = [4., 4., 4.]
    nst = len(heights)
    n_bays = len(spans_x)

    # Material properties
    eps_y = fy / elastic_modulus_steel
    Ec = (3320 * fc**0.5 + 6900) * 1000.
    n_seismic = None
    n_gravity = None
    configuration = "space"
    site = None
    site_idx = 0

    inputs = {"loads": [], "seismic": []}

    masses = None

    def __init__(self, flag3d=False):
        """
        initializes the input functions
        """
        self.flag3d = flag3d

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
        self.masses = np.zeros(len(self.heights))
        for st in range(self.nst):
            self.masses[st] = area / 9.81 * self.inputs["seismic"][st]

        return self.masses
