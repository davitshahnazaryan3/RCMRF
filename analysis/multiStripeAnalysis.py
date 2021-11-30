import multiprocessing as mp

import openseespy.opensees as op
import numpy as np
import pandas as pd
import os
from analysis.solutionAlgorithm import SolutionAlgorithm
from analysis.static import Static
from client.model import Model
from utils.utils import read_text_file, createFolder, export_to


def get_ground_motion(path):
    names_x = list(pd.read_csv(path / "GMR_H1_names.txt", header=None)[0])
    names_y = list(pd.read_csv(path / "GMR_H2_names.txt", header=None)[0])
    dts_list = read_text_file(path / "GMR_dts.txt")
    return names_x, names_y, dts_list


def get_ground_motion_batches(gm_dir):
    return next(os.walk(gm_dir))[1]


def get_records(gm_dir, output_dir):
    """
    Gets records for multiprocessing
    :param gm_dir: str
    :param output_dir: str          Output directory
    :return: None
    """
    # Get the ground motion set information
    gm_paths = get_ground_motion_batches(gm_dir)

    # Records in batches for multiprocessing
    records = {}

    # Set up directories to export MSA outputs for each batch of records
    for path in gm_paths:
        createFolder(output_dir / path)

        # Get the ground motion information
        names_x, names_y, dts_list = get_ground_motion(gm_dir / path)
        records[path] = {"X": names_x, "Y": names_y, "dt": dts_list}

    return records


class MultiStripeAnalysis:
    def __init__(self, sections_file, loads_file, materials, gm_dir, damping, omegas, output_path,
                 analysis_time_step=0.01, drift_capacity=10., system="perimeter", hingeModel="haselton", flag3d=False,
                 export_at_each_step=True, pflag=True):
        self.sections_file = sections_file
        self.loads_file = loads_file
        self.materials = materials
        self.gm_dir = gm_dir
        self.damping = damping
        self.omegas = omegas
        self.output_path = output_path
        self.analysis_time_step = analysis_time_step
        self.drift_capacity = drift_capacity
        self.system = system
        self.hingeModel = hingeModel.lower()
        self.flag3d = flag3d
        self.export_at_each_step = export_at_each_step
        self.pflag = pflag

        # Constants
        self.g = 9.81
        self.EXTRA_DUR = 10.0
        self.PTAGX = 10
        self.PTAGY = 20
        self.TSTAGX = 51
        self.TSTAGY = 52

        # Storing results
        self.outputs = {}
        self.IM_output = {}

    def call_model(self, generate_model=True):
        """
        Calls the Model
        Generates the model, defines gravity loads and defines time series
        :param generate_model: bool                         Generate model or not
        :return: class                                      Object Model
        """
        m = Model("MSA", self.sections_file, self.loads_file, self.materials, None, self.system,
                  self.hingeModel, flag3d=self.flag3d)
        if generate_model:
            # Create the nonlinear model
            m.model()
            # Define gravity loads
            m.define_loads(m.elements, apply_point=False)
            # Run static analysis
            s = Static()
            s.static_analysis(None, self.flag3d)
        return m

    def time_series(self, dt, pathx, pathy, fx, fy):
        """
        Defines time series
        :param dt: float                            Time step of record
        :param pathx: str                           Path to record of x direction
        :param pathy: str                           Path to record of y direction
        :param fx: float                            Scaling factor in x direction
        :param fy: float                            Scaling factor in y direction
        :return:
        """
        # Delete the old analysis and all it's component objects
        op.wipeAnalysis()

        w1 = self.omegas[0]
        w2 = self.omegas[1]
        a0 = 2 * w1 * w2 / (w2 ** 2 - w1 ** 2) * (w2 * self.damping - w1 * self.damping)
        a1 = 2 * w1 * w2 / (w2 ** 2 - w1 ** 2) * (-self.damping / w2 + self.damping / w1)

        # Rayleigh damping
        op.rayleigh(a0, 0.0, 0.0, a1)
        op.timeSeries('Path', self.TSTAGX, '-dt', dt, '-filePath', str(pathx), '-factor', fx)
        op.timeSeries('Path', self.TSTAGY, '-dt', dt, '-filePath', str(pathy), '-factor', fy)
        op.pattern('UniformExcitation', self.PTAGX, 1, '-accel', self.TSTAGX)
        if self.flag3d:
            op.pattern('UniformExcitation', self.PTAGY, 2, '-accel', self.TSTAGY)

        if self.flag3d:
            op.constraints('Penalty', 1.0e15, 1.0e15)
            # op.constraints("Transformation")
        else:
            op.constraints('Plain')
        op.numberer('RCM')
        op.system('UmfPack')

    def start_process(self, records):
        # Get number of CPUs available
        workers = mp.cpu_count()

        with mp.Pool(workers - 1) as pool:
            outputs = pool.imap(self.run_msa, list(records.items()))

            for _ in outputs:
                print("[SUCCESS]")

    def run_msa(self, item):
        name, data = item
        print(f"[START] Running {name} records...")

        names_x = data["X"]
        names_y = data["Y"]
        dts = data["dt"]

        # Initialize outputs
        self.outputs[name] = {}

        # For each record pair
        for rec in range(len(names_x)):
            eq_name_x = self.gm_dir / name / names_x[rec]
            eq_name_y = self.gm_dir / name / names_y[rec]
            dt = dts[rec]

            # duration
            accg_x = read_text_file(eq_name_x)
            dur = self.EXTRA_DUR + dt * len(accg_x)

            # Create the model
            m = self.call_model()

            # Create the time series
            self.time_series(dt, eq_name_x, eq_name_y, 1, 1)

            if self.pflag:
                print(f"[MSA] Record: {rec} - {name}: {names_x[rec]} and {names_y[rec]} pair;")

            th = SolutionAlgorithm(self.analysis_time_step, dur, self.drift_capacity, m.g.tnode, m.g.bnode, self.pflag,
                                   self.flag3d)
            self.outputs[name][rec] = th.ntha_results

            if self.export_at_each_step:
                export_to(self.output_path / "MSA" / name / f"Record{rec + 1}", self.outputs[name][rec], "pickle")

            # Wipe the model
            op.wipe()

        return self.outputs
