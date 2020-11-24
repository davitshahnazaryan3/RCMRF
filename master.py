"""
Master file to run OpenSees
"""
import openseespy.opensees as ops
from client.model import Model
from pathlib import Path
import os
import json
import pickle
from analysis.ida_htf import IDA_HTF


class Master:
    def __init__(self, sections_file, loads_file, materials_file, outputsDir, gmdir=None, gmfileNames=None, IM_type=2,
                 max_runs=15, analysis_time_step=.01, drift_capacity=20, analysis_type=None, system="Perimeter",
                 hingeModel="Hysteretic"):
        """
        Initializes master file
        :param sections_file: str                   Name of file containing section data in '*.csv' format
        :param loads_file: str                      Name of file containing load data in '*.csv' format
        :param materials_file: dict                 Concrete and reinforcement material properties
        :param outputsDir: str
        :param gmdir: str
        :param gmfileNames: list(str)
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
        :param hingeModel: str                      Hinge model (hysteretic or haselton)
        """
        self.sections_file = sections_file
        self.loads_file = loads_file
        self.materials_file = materials_file
        self.outputsDir = outputsDir
        self.gmdir = gmdir
        self.gmfileNames = gmfileNames
        self.IM_type = IM_type
        self.max_runs = max_runs
        self.analysis_time_step = analysis_time_step
        self.drift_capacity = drift_capacity
        self.analysis_type = analysis_type
        self.system = system
        self.APPLY_GRAVITY_ELF = False
        self.FIRST_INT = .05
        self.INCR_STEP = .05
        self.DAMPING = .05
        # For IDA
        self.NAME_X_FILE = gmdir / gmfileNames[0]
        self.NAME_Y_FILE = gmdir / gmfileNames[1]
        self.DTS_FILE = gmdir / gmfileNames[2]
        self.DURS_FILE = gmdir / gmfileNames[3]
        self.hingeModel = hingeModel
        # Create an outputs directory if none exists
        self.createFolder(outputsDir)

    def createFolder(self, directory):
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
        ops.wipe()

    def call_model(self, generate_model=True):
        """
        Calls the Model
        :param generate_model: bool                         Generate model or not
        :return: class                                      Object Model
        """
        m = Model(self.analysis_type, self.sections_file, self.loads_file, self.materials_file, self.outputsDir,
                  self.system, hingeModel=self.hingeModel)
        if generate_model:
            m.model()
        return m

    def run_model(self):
        """
        Initializes model creator and runs analysis
        :return: None
        """
        # Create model
        if "TH" in self.analysis_type or "timehistory" in self.analysis_type or "IDA" in self.analysis_type:
            m = self.call_model(generate_model=False)
        else:
            m = self.call_model()

        # Run analysis and record the goodies
        if self.analysis_type is None or -1 in self.analysis_type or not self.analysis_type:
            self.analysis_type = []
            print("[CHECK] Checking the integrity of the model")
            m.perform_analysis()
            print("[SUCCESS] Integrity of model verified")

        elif "PO" in self.analysis_type or "pushover" in self.analysis_type:

            try:
                with open(self.outputsDir / "MA.json") as f:
                    modal_analysis_outputs = json.load(f)
                m.perform_analysis(spo_pattern=2, mode_shape=modal_analysis_outputs["Mode1"])

            except:
                print("[WARNING] 1st Mode-proportional loading was selected. MA outputs are missing! "
                      "\nRunning with triangular pattern")

                m.perform_analysis(spo_pattern=1)

        elif "ELF" in self.analysis_type or "ELFM" in self.analysis_type and self.APPLY_GRAVITY_ELF:
            m.perform_analysis()

        elif "TH" in self.analysis_type or "timehistory" in self.analysis_type or "IDA" in self.analysis_type:
            print("[INITIATE] IDA started")
            try:
                with open(self.outputsDir / "MA.json") as f:
                    results = json.load(f)
                    period = results["Periods"][0]
                    damping = results["Damping"][0]
                    omegas = results["CircFreq"]
            except:
                raise ValueError("[EXCEPTION] Static and modal analysis data do not exist.")

            # Initialize
            ida = IDA_HTF(self.FIRST_INT, self.INCR_STEP, self.max_runs, self.IM_type, period, damping, omegas,
                          self.analysis_time_step, self.drift_capacity, self.NAME_X_FILE, self.NAME_Y_FILE,
                          self.DTS_FILE, self.DURS_FILE, self.gmdir, self.analysis_type, self.sections_file,
                          self.loads_file, self.materials_file, self.system, hingeModel=self.hingeModel,
                          pflag=True)

            # Set-up
            ida.establish_im(m.g.tnode, m.g.bnode)

            # Export results
            with open(self.outputsDir / "IDA.pickle", "wb") as handle:
                pickle.dump(ida.outputs, handle)

            print("[SUCCESS] IDA done")

        else:
            m.define_loads(m.elements)
            m.perform_analysis(damping=self.DAMPING)

        # Record outputs
        if not bool(m.results):
            pass
        else:
            directory = self.outputsDir / "results.pkl"
            with open(directory, "wb") as handle:
                pickle.dump(m.results, handle)


if __name__ == "__main__":

    def truncate(n, decimals=0):
        multiplier = 10 ** decimals
        return int(n * multiplier) / multiplier

    def get_time(start_time):
        elapsed = timeit.default_timer() - start_time
        print('Running time: ', truncate(elapsed, 2), ' seconds')
        print('Running time: ', truncate(elapsed / float(60), 2), ' minutes')

    import timeit
    directory = Path.cwd().parents[0] / ".applications/case1/Output"

    start_time = timeit.default_timer()

    materials_file = directory / "materials.csv"
    section_file = directory / "hinge_models.csv"
    loads_file = directory / "action.csv"
    analysis_type = ["ST"]
    hingeModel = "Hysteretic"
    outputsDir = directory / "RCMRF/archive"
    # GM directory
    gmdir = Path.cwd() / "groundMotion"
    gmfileNames = ["GMR_names1.txt", "GMR_names2.txt", "GMR_dts.txt", "GMR_durs.txt"]

    m = Master(section_file, loads_file, materials_file, outputsDir, gmdir=gmdir, gmfileNames=gmfileNames,
               analysis_type=analysis_type, system="Perimeter",  hingeModel=hingeModel)
    m.wipe()
    m.run_model()
    m.wipe()

    # get_time(start_time)
