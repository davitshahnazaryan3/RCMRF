"""
Master file to run OpenSees
3D model may be used with hingeModel = Hysteretic only - haselton to be adapted
It is recommended to run each analysis separately with the order being:
1. Static analysis - ST
2. Modal analysis - MA
3. Static pushover analysis - PO
4. Nonlinear time history analysis - TH (supports IDA only for now)

Requires 3-4 input files
Materials file for concrete and reinforcement properties
Action file defining gravity loads and pdelta loads (not necessary for 3D)
Hinge model information in primary direction
Hinge model information in secondary direction (not necessary for 2D)
Hinge model information for gravity frame elements (not necessary for 2D)
"""
import openseespy.opensees as op
from client.model import Model
from pathlib import Path
import numpy as np
import os
import json
import pickle
from analysis.ida_htf_3d import IDA_HTF_3D


class Master:
    def __init__(self, sections_file, loads_file, materials_file, outputsDir, gmdir=None, gmfileNames=None, IM_type=2,
                 max_runs=15, analysis_time_step=.01, drift_capacity=10., analysis_type=None, system="Perimeter",
                 hinge_model="Hysteretic", flag3d=False, direction=0, export_at_each_step=False,
                 period_assignment=None, periods_ida=None):
        """
        Initializes master file
        :param sections_file: str                   Name of file containing section data in '*.csv' format
        :param loads_file: str                      Name of file containing load data in '*.csv' format
        :param materials_file: dict                 Concrete and reinforcement material properties
        :param outputsDir: str                      Directory for exporting the analyses results
        :param gmdir: str                           Ground motion directory
        :param gmfileNames: list(str)               Files containing information on ground motion files (check sample)
        :param IM_type: int                         Intensity measure tag for IDA
        :param max_runs: int                        Maximum number of runs for IDA
        :param analysis_time_step: float            Analysis time step for IDA
        :param drift_capacity: float                Drift capacity in % for stopping IDA
        :param analysis_type: list(str)             Type of analysis for which we are recording [TH, PO, ST, MA, ELF]
                                                    TH - time history analysis
                                                    PO - static pushover analysis
                                                    ST - static analysis (e.g. gravity)
                                                    MA - modal analysis
                                                    ELF - equivalent lateral force method of analysis
        :param system: str                          System type (e.g. Perimeter or Space)
        :param hinge_model: str                     Hinge model (hysteretic or haselton)
        :param flag3d: bool                         True for 3D modelling, False for 2D modelling
        :param direction: int                       Direction of analysis
        :param export_at_each_step: bool            Exporting at each step, i.e. record-run
        :param period_assignment: dict              Period assigment IDs (for 3D only)
        :param periods_ida: list                    Periods to use for IDA, optional, in case MA periods are not needed
        """
        # TODO, add support for Haselton, currently only a placeholder, need to adapt for 3D etc.
        # list of strings for 3D modelling, and string for 2D modelling
        self.sections_file = sections_file

        # Input arguments
        self.outputsDir = outputsDir
        self.loads_file = loads_file
        self.materials_file = materials_file
        self.gmdir = gmdir
        self.gmfileNames = gmfileNames
        self.IM_type = IM_type
        self.max_runs = max_runs
        self.analysis_time_step = analysis_time_step
        self.drift_capacity = drift_capacity
        self.analysis_type = analysis_type
        self.system = system
        self.hinge_model = hinge_model.lower()
        self.flag3d = flag3d
        self.export_at_each_step = export_at_each_step
        self.direction = direction
        self.period_assignment = period_assignment
        self.periods_ida = periods_ida
        self.APPLY_GRAVITY_ELF = False
        self.FIRST_INT = .05
        self.INCR_STEP = .05
        self.DAMPING = .05

        # For IDA
        self.NAME_X_FILE = gmdir / gmfileNames[0]
        self.NAME_Y_FILE = gmdir / gmfileNames[1]
        self.DTS_FILE = gmdir / gmfileNames[2]
        self.DURS_FILE = gmdir / gmfileNames[3]

        # Create an outputs directory if none exists
        self.createFolder(outputsDir)

        if direction != 0 and flag3d and analysis_type == "MA":
            print("[WARNING] Direction should be set to 0 for Modal Analysis!")

    @staticmethod
    def createFolder(directory):
        """
        Checks whether provided directory exists, if no creates one
        :param directory: str
        :return: None
        """
        try:
            if not os.path.exists(directory):
                os.makedirs(directory)
        except OSError:
            print('Error: Creating directory. ' + directory)

    @staticmethod
    def wipe():
        """
        Perform a clean wipe
        :return: None
        """
        op.wipe()

    def call_model(self, generate_model=True):
        """
        Calls the Model
        :param generate_model: bool                         Generate model or not
        :return: class                                      Object Model
        """
        m = Model(self.analysis_type, self.sections_file, self.loads_file, self.materials_file, self.outputsDir,
                  self.system, hingeModel=self.hinge_model, flag3d=self.flag3d, direction=self.direction)
        # Generate the model if specified
        if generate_model:
            m.model()
        return m

    def run_model(self):
        """
        Initializes model creator and runs analysis
        :return: None
        """
        # Run analysis and record the goodies
        if self.analysis_type is None or -1 in self.analysis_type or not self.analysis_type:
            # Checks the integrity of the files for any sources of error
            self.analysis_type = []
            print("[CHECK] Checking the integrity of the model")
            m = self.call_model()
            m.perform_analysis()
            print("[SUCCESS] Integrity of model verified")

        elif "PO" in self.analysis_type or "pushover" in self.analysis_type:
            try:
                # Static pushover analysis needs to be run after modal analysis
                with open(self.outputsDir / "MA.json") as f:
                    modal_analysis_outputs = json.load(f)
                # Modal shape as the SPO lateral load pattern shape
                if self.direction == 0:
                    mode_shape = np.abs(modal_analysis_outputs["Mode1"])
                else:
                    mode_shape = np.abs(modal_analysis_outputs["Mode2"])
                # Normalize, helps to avoid convergence issues
                mode_shape = mode_shape / max(mode_shape)
                mode_shape = np.round(mode_shape, 2)
                # Call and run the OpenSees model
                m = self.call_model()
                m.perform_analysis(spo_pattern=2, mode_shape=mode_shape)

            except:
                print("[WARNING] 1st Mode-proportional loading was selected. MA outputs are missing! "
                      "\nRunning with triangular pattern")
                m = self.call_model()
                m.perform_analysis(spo_pattern=1)

        elif "ELF" in self.analysis_type or "ELFM" in self.analysis_type and self.APPLY_GRAVITY_ELF:
            if self.flag3d:
                # It might still run though :)
                raise ValueError("[EXCEPTION] ELF for a 3D model not supported!")

            # Equivalent lateral force method of analysis
            m = self.call_model()
            m.perform_analysis()

        elif "TH" in self.analysis_type or "timehistory" in self.analysis_type or "IDA" in self.analysis_type:
            # Create a folder for NLTHA
            self.createFolder(self.outputsDir / "NLTHA")

            # Nonlinear time history analysis
            print("[INITIATE] IDA started")
            try:
                if self.periods_ida is None:
                    with open(self.outputsDir / "MA.json") as f:
                        results = json.load(f)
                    if self.flag3d:
                        period = [results["Periods"][self.period_assignment["x"]],
                                  results["Periods"][self.period_assignment["y"]]]
                    else:
                        period = results["Periods"][0]
                    damping = results["Damping"][0]
                    omegas = results["CircFreq"]
                else:
                    period = self.periods_ida
                    damping = 0.05
                    omegas = 2 * np.pi / (np.array(period))

            except:
                raise ValueError("[EXCEPTION] Modal analysis data does not exist.")

            # Initialize
            ida = IDA_HTF_3D(self.FIRST_INT, self.INCR_STEP, self.max_runs, self.IM_type, period, damping, omegas,
                             self.analysis_time_step, self.drift_capacity, self.NAME_X_FILE, self.NAME_Y_FILE,
                             self.DTS_FILE, self.DURS_FILE, self.gmdir, self.analysis_type, self.sections_file,
                             self.loads_file, self.materials_file, self.system, hingeModel=self.hinge_model,
                             pflag=True, flag3d=self.flag3d, export_at_each_step=self.export_at_each_step)

            # The Set-up
            ida.establish_im(output_dir=self.outputsDir / "NLTHA")

            # Export results
            if self.export_at_each_step:
                with open(self.outputsDir / "IDA.pickle", "wb") as handle:
                    pickle.dump(ida.outputs, handle)
            if os.path.exists(self.outputsDir / "IM.csv"):
                im_filename = "IM_temp.csv"
            else:
                im_filename = "IM.csv"
            np.savetxt(self.outputsDir / im_filename, ida.IM_output, delimiter=',')

            print("[SUCCESS] IDA done")

        else:
            # Runs static or modal analysis (ST or MA)
            m = self.call_model()
            m.define_loads(m.elements, apply_point=False)
            m.perform_analysis(damping=self.DAMPING)


if __name__ == "__main__":

    def truncate(n, decimals=0):
        multiplier = 10 ** decimals
        return int(n * multiplier) / multiplier

    def get_time(start_time):
        elapsed = timeit.default_timer() - start_time
        print('Running time: ', truncate(elapsed, 2), ' seconds')
        print('Running time: ', truncate(elapsed / float(60), 2), ' minutes')

    import timeit
    start_time = timeit.default_timer()

    # Directories
    directory_x = Path.cwd().parents[0] / ".applications/LOSS Validation Manuscript/Case21/Cache/framex"
    materials_file = directory_x.parents[1] / "materials.csv"
    loads_file = directory_x.parents[1] / "action.csv"
    outputsDir = directory_x.parents[1] / "RCMRF"

    # X direction
    section_file_x = directory_x / "hinge_models.csv"

    # Y direction
    directory_y = Path.cwd().parents[0] / ".applications/LOSS Validation Manuscript/Case21/Cache/framey"
    section_file_y = directory_y / "hinge_models.csv"

    # Gravity elements
    directory_gr = Path.cwd().parents[0] / ".applications/LOSS Validation Manuscript/Case21/Cache"
    section_file_gr = directory_gr / "gravity_hinges.csv"

    # Directories
    section_file = {"x": section_file_x, "y": section_file_y, "gravity": section_file_gr}

    # GM directory
    gmdir = Path.cwd() / "sample/groundMotion"
    gmfileNames = ["GMR_names1.txt", "GMR_names2.txt", "GMR_dts.txt", "GMR_durs.txt"]

    # RCMRF inputs
    hingeModel = "Hysteretic"
    analysis_type = ["TH"]
    # Run MA always with direction = 0
    direction = 0
    flag3d = True
    export_at_each_step = True
    period_assignment = {"x": 0, "y": 1}
    periods = [0.72, 0.62]
    # Let's go...
    m = Master(section_file, loads_file, materials_file, outputsDir, gmdir=gmdir, gmfileNames=gmfileNames,
               analysis_type=analysis_type, system="Perimeter", hinge_model=hingeModel, flag3d=flag3d,
               direction=direction, export_at_each_step=export_at_each_step, period_assignment=period_assignment,
               periods_ida=periods, max_runs=15)

    m.wipe()
    m.run_model()

    # Wipe the model
    m.wipe()

    # Time it
    get_time(start_time)
