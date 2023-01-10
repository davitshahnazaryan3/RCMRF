import multiprocessing as mp
import openseespy.opensees as op
import pandas as pd
import numpy as np
import os

from analysis.solutionAlgorithm import SolutionAlgorithm
from analysis.static import Static
from client.model import Model
from utils.utils import read_text_file, createFolder, export_to


def get_ground_motion(path, gmfileNames):
    names_y = None

    if len(gmfileNames) == 2:
        # 2D modelling, a single direction of a ground motion record is provided
        names_x = list(pd.read_csv(path / gmfileNames[0], header=None)[0])
        dts_list = read_text_file(path / gmfileNames[1])
    else:
        names_x = list(pd.read_csv(path / gmfileNames[0], header=None)[0])
        names_y = list(pd.read_csv(path / gmfileNames[1], header=None)[0])
        dts_list = read_text_file(path / gmfileNames[2])

    try:
        dts_list = [float(dts_list)]
    except:
        dts_list = list(dts_list)

    return names_x, names_y, dts_list


def get_ground_motion_batches(gm_dir):
    return next(os.walk(gm_dir))[1]


def get_records(gm_dir, output_dir, gmfileNames):
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
        names_x, names_y, dts_list = get_ground_motion(gm_dir / path, gmfileNames)
        records[path] = {"X": names_x, "Y": names_y, "dt": dts_list}

    return records


class MultiStripeAnalysis:
    use_multiprocess = False

    def __init__(self, sections_file, loads_file, materials, gm_dir, damping, omegas, output_path,
                 analysis_time_step=None, drift_capacity=10., system="space", hingeModel="hysteretic", flag3d=False,
                 export_at_each_step=True, pflag=True, gm_scaling_factor=1., use_recorder=True, recorder_cache=None):
        """
        Initialize
        :param sections_file: Path or dict
        :param loads_file: Path
        :param materials: Path
        :param gm_dir: Path
        :param damping: float
        :param omegas: List(float)
        :param output_path: Path
        :param analysis_time_step: float
        :param drift_capacity: float
        :param system: str
        :param hingeModel: str
        :param flag3d: bool
        :param export_at_each_step: bool
        :param pflag: bool
        :param gm_scaling_factor: float
        :param use_recorder: bool                   Uses openseespy recorder to output file instead of node recorders
        :param recorder_cache: str
        """
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
        self.gm_scaling_factor = gm_scaling_factor
        self.use_recorder = use_recorder
        self.recorder_cache = recorder_cache

        # Constants
        self.g = 9.81
        self.EXTRA_DUR = 10.0
        self.PTAGX = 10
        self.PTAGY = 20
        self.TSTAGX = 51
        self.TSTAGY = 52

        # Storing results
        self.outputs = {}

    def _call_model(self, generate_model=True):
        """
        Calls the Model
        Generates the model, defines gravity loads and defines time series
        :param generate_model: bool                         Generate model or not
        :return: class                                      Object Model
        """
        m = Model("MSA", self.sections_file, self.loads_file, self.materials, self.output_path, self.system,
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

    @staticmethod
    def _append_record(x, y):
        """
        Record selection issue related to the database of records
        Whereby the two pairs of the records have different sizes
        To remedy the issue, the smallest of the records is appended with zeros
        :param x: list
        :param y: list
        :return: list, list
        """
        nx = len(x)
        ny = len(y)

        if nx < ny:
            x = np.append(x, np.zeros(ny - nx))
        if ny < nx:
            y = np.append(y, np.zeros(nx - ny))
        return x, y

    def _time_series(self, dt, pathx, pathy, fx, fy):
        """
        Defines time series
        :param dt: float                            Time step of record
        :param pathx: str                           Path to record of x direction
        :param pathy: str                           Path to record of y direction
        :param fx: float                            Scaling factor in x direction
        :param fy: float                            Scaling factor in y direction
        :return: None
        """
        # Eigen value analysis
        # eigen_values = np.asarray(op.eigen(max(self.damping)))
        # omega = eigen_values ** 0.5
        # periods = 2 * np.pi / omega
        # self.omegas = omega

        # Delete the old analysis and all it's component objects
        op.wipeAnalysis()

        w1 = self.omegas[0]
        w2 = self.omegas[1]
        a0 = 2 * w1 * w2 / (w2 ** 2 - w1 ** 2) * (w2 * self.damping - w1 * self.damping)
        a1 = 2 * w1 * w2 / (w2 ** 2 - w1 ** 2) * (-self.damping / w2 + self.damping / w1)

        # Rayleigh damping
        op.rayleigh(a0, 0.0, 0.0, a1)
        op.timeSeries('Path', self.TSTAGX, '-dt', dt, '-filePath', str(pathx), '-factor', fx)
        if pathy is not None:
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

    def start_process(self, records, workers=0):
        """
        Start the parallel computation
        :param records: dict
        :param workers: int
        :return: None
        """
        self.use_multiprocess = True

        # Get number of CPUs available
        if workers == 0:
            workers = mp.cpu_count()
        if workers > 0:
            workers = workers + 1

        with mp.Pool(workers - 1, maxtasksperchild=1) as pool:
            outputs = pool.imap(self.run_msa, list(records.items()))

            for _ in outputs:
                print("[SUCCESS]")

    def run_msa(self, item):
        """
        Records
        :param item: List(str, dict)
        :return: dict
        """
        name, data = item
        print(f"[START] Running {name} records...")

        # Get the filenames of record pairs and time steps
        names_x = data["X"]
        names_y = data["Y"]
        dts = data["dt"]

        # initialize for second direction (if 3D modelling is utilized, both directions of record, then the variables
        # will be updated)
        accg_y = None
        eq_name_y = None

        # Initialize outputs
        self.outputs[name] = {}

        # For each record pair
        for rec in range(len(names_x)):
            if self.use_multiprocess:
                self.recorder_cache = f"{name}_{rec}.txt"

            # reading records
            eq_name_x = self.gm_dir / name / names_x[rec]
            dt = dts[rec]
            accg_x = read_text_file(eq_name_x)

            # Second direction
            if names_y is not None:
                eq_name_y = self.gm_dir / name / names_y[rec]
                accg_y = read_text_file(eq_name_y)
                # duration, make sure both directions have the same size
                accg_x, accg_y = self._append_record(accg_x, accg_y)

            # add extra duration of free vibrations to the records
            dur = round(self.EXTRA_DUR + dt * len(accg_x), 5)

            # analysis time step
            if self.analysis_time_step is None:
                analysis_time_step = dt
            else:
                analysis_time_step = self.analysis_time_step

            # Create the model
            m = self._call_model()

            # Create the time series
            self._time_series(dt, eq_name_x, eq_name_y, self.gm_scaling_factor, self.gm_scaling_factor)

            if self.pflag:
                if names_y is None:
                    print(f"[MSA] Record: {rec} - {name}: {names_x[rec]};")
                else:
                    print(f"[MSA] Record: {rec} - {name}: {names_x[rec]} and {names_y[rec]} pair;")

            analysis_time_step = min(analysis_time_step, dt)
            if dt % analysis_time_step != 0:
                analysis_time_step = dt / (int(dt / analysis_time_step))

            th = SolutionAlgorithm(analysis_time_step, dur, self.drift_capacity, m.g.tnode, m.g.bnode,
                                   dt, accg_x, accg_y, self.gm_scaling_factor, self.gm_scaling_factor,
                                   pflag=self.pflag, flag3d=self.flag3d, use_recorder=self.use_recorder,
                                   recorder_cache=self.recorder_cache)
            self.outputs[name][rec] = th.ntha_results

            if self.export_at_each_step:
                export_to(self.output_path / "MSA" / name / f"Record{rec + 1}", self.outputs[name][rec], "pickle")

            # Wipe the model
            op.wipe()

        return self.outputs
