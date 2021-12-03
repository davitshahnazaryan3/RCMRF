"""
Incremental Dynamic Analysis using Hunt, Trace and Fill algorithm
Object temporarily kept. Will be deleted later...
"""
import openseespy.opensees as op
import numpy as np
import pandas as pd
from analysis.solutionAlgorithm import SolutionAlgorithm
from analysis.static import Static
from client.model import Model


class IDA_HTF:
    def __init__(self, first_int, incr_step, max_runs, IM_type, T_info, xi, omegas, dt, dcap, nmsfile_x, nmsfile_y,
                 dts_file, durs_file, gm_dir, analysis_type, sections_file, loads_file, materials, system='space',
                 hingeModel='hysteretic', pflag=True, flag3d=False):
        """
        Initializes IDA
        :param first_int: float                     The first intensity to run the elastic run (e.g. 0.05g)
        :param incr_step: float                     The increment used during the hunting phase (e.g. 0.10g)
        :param max_runs: int                        Maximum number of runs to use (e.g. 20)
        :param IM_type: int                         Intensity measure with which to conduct IDA
                                                    1. PGA
                                                    2. Sa at a given period (e.g. Sa(T1)), geometric mean of two record
                                                    Sa(T) as the IM
                                                    3. more to be added
        :param T_info: list                         List of period info required by specified IM
                                                    1. [], will ignore any entries if present
                                                    2. [1.0], single value of period to condition to
                                                    3. more to be added
        :param xi: float                            Elastic damping, typically 0.05
        :param omegas: list                         Circular frequences
        :param dt: float                            Analysis time step
        :param dcap: float                          Drift capacity in %
        :param nmsfile_x: str                       Text file with names of the X direction records in the form "*.txt"
        :param nmsfile_y: str                       Text file with names of the Y direction records in the form "*.txt"
        :param dts_file: str                        Text file with the time steps of the records
        :param durs_file: str                       Text file with the durations of the records
        :param gm_dir: str                          Directory containg the Ground Motions
        :param analysis_type, sections_file, loads_file, materials, system, hingeModel: See client\model.py
        :param pflag: bool                          Whether print information on screen or not
        :param flag3d: bool                         True for 3D modelling, False for 2D modelling
        """
        self.first_int = first_int
        self.incr_step = incr_step
        self.max_runs = max_runs
        self.IM_type = IM_type
        self.T_info = T_info
        self.xi = xi
        self.omegas = omegas
        self.dt = dt
        self.dcap = dcap
        self.nmsfile_x = nmsfile_x
        self.nmsfile_y = nmsfile_y
        self.dts_file = dts_file
        self.durs_file = durs_file
        self.gm_dir = gm_dir
        self.g = 9.81
        self.EXTRA_DUR = 10.0
        self.analysis_type = analysis_type
        self.sections_file = sections_file
        self.loads_file = loads_file
        self.materials = materials
        self.hingeModel = hingeModel
        self.system = system
        self.pflag = pflag
        self.flag3d = flag3d
        self.PTAGX = 10
        self.PTAGY = 20
        self.TSTAGX = 51
        self.TSTAGY = 52
        self.outputs = {}
        # Initialize intensity measures
        self.IM_output = None

    def call_model(self, generate_model=True):
        """
        Calls the Model
        Generates the model, defines gravity loads and defines time series
        :param generate_model: bool                         Generate model or not
        :return: class                                      Object Model
        """
        m = Model(self.analysis_type, self.sections_file, self.loads_file, self.materials, None, self.system,
                  self.hingeModel, flag3d=self.flag3d)
        if generate_model:
            # Create the nonlinear model
            m.model()
            # Define gravity loads
            m.define_loads(m.elements, apply_point=True)
            # Run static analysis
            s = Static()
            s.static_analysis(self.flag3d)
        return m

    def read_text_file(self, name, col):
        """
        Reads a text file
        :param name: str                            Name of text file in "*.txt" format
        :param col: int                             Column of text file to read from
        :return: ndarray                            Read column in np.array format
        """
        f = open(name, 'r')
        lines = f.readlines()
        data = np.array([])
        for x in lines:
            data = np.append(data, float(x.split()[col]))
        f.close()
        return data

    def get_IM(self, eq, dt, T, xi):
        """
        Gets Sa(T) of a given record, for a specified value of period T using the Newmark Average Acceleration
        :param eq: str                              Filename which is a single column file in units of g (e.g. "Eq.txt")
        :param dt: float                            Time step in seconds (e.g. 0.01)
        :param T: float                             Period in seconds (e.g. 1.0)
        :param xi: float                            Elastic damping (e.g. 0.05)
        :return: floats                             sa - pseudo-spectral acceleration in g
                                                    sv - pseudo-spectral velocity in m/s
                                                    sd - spectral displacement in m
                                                    pga - peak ground acceleration in g
        """
        # Read the acceleration time series
        accg = self.read_text_file(eq, 0)

        if T == 0.0:
            pga = 0.0
            for i in range(len(accg)):
                temp2 = abs(accg[i])
                if temp2 > pga:
                    pga = temp2
            return pga

        else:
            # Initialize
            # Newmark terms (Set for Average Acceleration method)
            gamma = 0.5
            beta = 0.25
            # Set the mass to 1kg
            ms = 1.0

            # Change the units of the record to m/s2
            acc = accg * self.g
            # Create the force in N
            p = -ms * acc

            # Create time vector
            t = np.arange(0., (len(acc)-1)*dt, dt)

            # Calculate the initial values
            # Stiffness in N/m
            k = ms * np.power(6.2832 / T, 2)
            # Circular frequency
            w = np.power(k / ms, 0.5)
            # Damping coefficient
            c = 2 * xi * ms * w
            # Initial acceleration in m/s2
            a0 = p[0] / ms
            # Stiffness term (see Chopra book)
            k_bar = k + gamma * c / beta / dt + ms / beta / dt ** 2
            # Constant term (see Chopra book)
            A = ms / beta / dt + gamma * c / beta
            # Constant term (see Chopra book)
            B = ms / 2 / beta + dt * c * (gamma / 2 / beta - 1)

            # Initialize some vectors
            u = np.array([0.0])
            v = np.array([0.0])
            a = np.array([a0])
            du = np.array([])
            dv = np.array([])
            da = np.array([])
            dp = np.array([])
            dp_bar = np.array([])

            for i in range(len(t) - 1):
                dp = np.append(dp, p[i + 1] - p[i])
                dp_bar = np.append(dp_bar, dp[i] + A * v[i] + B * a[i])
                du = np.append(du, dp_bar[i] / k_bar)
                dv = np.append(dv, gamma * du[i] / beta / dt - gamma * v[i] / beta + dt * (1 - gamma / 2 / beta) * a[i])
                da = np.append(da, du[i] / beta / dt ** 2 - v[i] / beta / dt - a[i] / 2 / beta)
                u = np.append(u, u[i] + du[i])
                v = np.append(v, v[i] + dv[i])
                a = np.append(a, a[i] + da[i])

            umax = 0.0
            for i in range(len(t)):
                # Absolute current displacement
                temp1 = abs(u[i])
                if temp1 > umax:
                    umax = temp1

            pga = 0.0
            for i in range(len(accg)):
                # Absolute current ground acceleration
                temp2 = abs(accg[i])
                if temp2 > pga:
                    pga = temp2

            # Calculate spectral values
            sd = umax
            sv = sd * w
            sa = sd * w ** 2 / self.g
            return sd, sv, sa

    def get_gm(self):
        """
        Gets ground motion information (i.e. names, time steps, durations)
        :return: lists                              List of names, time steps and durations of each record
        """
        eqnms_list_x = list(pd.read_csv(self.nmsfile_x, header=None)[0])
        eqnms_list_y = list(pd.read_csv(self.nmsfile_y, header=None)[0])
        dts_list = self.read_text_file(self.dts_file, 0)
        durs_list = self.read_text_file(self.durs_file, 0)
        return eqnms_list_x, eqnms_list_y, dts_list, durs_list

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
        w1 = self.omegas[0]
        w2 = self.omegas[1]
        a0 = 2 * w1 * w2 / (w2 ** 2 - w1 ** 2) * (w2 * self.xi - w1 * self.xi)
        a1 = 2 * w1 * w2 / (w2 ** 2 - w1 ** 2) * (-self.xi / w2 + self.xi / w1)

        # Rayleigh damping
        op.rayleigh(a0, 0.0, 0.0, a1)
        op.timeSeries('Path', self.TSTAGX, '-dt', dt, '-filePath', str(pathx), '-factor', fx)
        # op.timeSeries('Path', self.TSTAGY, '-dt', dt, '-filePath', str(pathy), '-factor', fy)
        op.pattern('UniformExcitation', self.PTAGX, 1, '-accel', self.TSTAGX)

        # Delete the old analysis and all it's component objects
        op.wipeAnalysis()
        op.constraints('Plain')
        op.numberer('Plain')
        op.system('UmfPack')

    def establish_im(self):
        """
        Establishes IM and performs analyses
        :return: None
        """
        # Get the ground motion set information
        eqnms_list_x, eqnms_list_y, dts_list, durs_list = self.get_gm()
        nrecs = len(dts_list)

        # Initialize intensity measures (shape)
        self.IM_output = np.zeros((nrecs, self.max_runs))

        # Loop for each record
        for rec in range(nrecs):
            # Counting starts from 0
            self.outputs[rec] = {}
            # Get the ground motion set information
            eq_name_x = self.gm_dir / eqnms_list_x[rec]
            eq_name_y = self.gm_dir / eqnms_list_y[rec]
            dt = dts_list[rec]
            dur = self.EXTRA_DUR + durs_list[rec]

            # Establish the IM
            if self.IM_type == 1:
                print('[IDA] IM is the PGA')
                pga = self.get_IM(eq_name_x, dts_list[rec], 0.0, self.xi)
                IMx = pga
                pga = self.get_IM(eq_name_y, dts_list[rec], 0.0, self.xi)
                IMy = pga
                # Get the geometric mean
                IM_geomean = np.power(IMx * IMy, 0.5)

            elif self.IM_type == 2:
                print('[IDA] IM is Sa at a specified period')
                Tcond = self.T_info
                sd, sv, sa = self.get_IM(eq_name_x, dts_list[rec], Tcond, self.xi)
                IMx = sa
                sd, sv, sa = self.get_IM(eq_name_y, dts_list[rec], Tcond, self.xi)
                IMy = sa
                # Get the geometric mean
                IM_geomean = np.power(IMx * IMy, 0.5)

            else:
                raise ValueError('[EXCEPTION] IM type provided incorrectly (must be 1 or 2)')

            # Set up the initial indices for HTF
            j = 1
            IM = np.zeros((self.max_runs, ))                # Initialize the list of IM used
            IMlist = np.array([])                           # A list to be used in filling
            hflag = 1                                       # Hunting flag (1 for when we are hunting)
            tflag = 0                                       # Tracing flag (0 at first)
            fflag = 0                                       # Filling flag (0 at first)
            jhunt = 0                                       # The run number we hunted to

            # The aim is to run NLTHA max_runs times
            while j <= self.max_runs:
                # As long as hunting flag is 1, meaning that collapse have not been reached
                if hflag == 1:
                    print("[STEP] We join the hunt...")
                    # Determine the intensity to run at during the hunting
                    if j == 1:
                        # First IM
                        IM[j - 1] = self.first_int
                    else:
                        # Subsequent IMs
                        IM[j - 1] = IM[j - 2] + (j - 1) * self.incr_step

                    # Determine the scale factor that needs to be applied to the record
                    sf_x = IM[j - 1] / IM_geomean * self.g
                    sf_y = IM[j - 1] / IM_geomean * self.g

                    run = f"Record{rec + 1}_Run{j}"
                    if self.pflag:
                        print(f"[IDA] IM = {IM[j - 1]}")

                    # The hunting intensity has been determined, now analysis commences
                    m = self.call_model()
                    self.time_series(dt, eq_name_x, eq_name_y, sf_x, sf_y)
                    if self.pflag:
                        print(f"[IDA] Record: {rec + 1}; Run: {j}; IM: {IM[j - 1]}")

                    # Commence analysis
                    th = SolutionAlgorithm(self.dt, dur, self.dcap, m.g.tnode, m.g.bnode, self.pflag, self.flag3d)
                    self.outputs[rec][j] = th.ntha_results

                    # Check the hunted run for collapse
                    """
                    c_index = -1                Analysis failed to converge at controltime of Tmax
                    c_index = 0                 Analysis completed successfully
                    c_index = 1                 Local structure collapse
                    """
                    # Check the hunted run for collapse
                    if th.c_index > 0:
                        # Collapse is caught, so stop hunting
                        hflag = 0
                        # Start tracing
                        tflag = 1
                        # The value of j we hunted to
                        jhunt = j
                        # Check whether first increment is too large
                        if jhunt == 2:
                            print(f"[IDA] Warning: {run} - Collapsed achieved on first increment, reduce increment...")
                    else:
                        self.IM_output[rec, j-1] = IM[j - 1]
                        # increment j
                        j += 1
                    op.wipe()

                # When first collapse is reached, tracing commences between last convergence and the first collapse
                if tflag == 1:
                    print("[STEP] Tracing...")
                    # First phase is to trace the last DeltaIM to get within the resolution
                    if j == jhunt:
                        # This is the IM of the hunting collapse (IM to cause collapse)
                        firstC = IM[j - 1]
                        # Remove that value of IM from the array
                        IM[j - 1] = 0.0

                    # Determine the difference between the hunting's noncollapse and collapse IM
                    diff = firstC - IM[j - 2]

                    # Take 0.2 of the difference
                    inctr = 0.2 * diff

                    # Place a lower threshold on the increment so it doesnt start tracing too fine
                    if inctr < 0.05:
                        inctr = 0.025

                    # Calculate new tracing IM, which is previous noncollapse plus increment
                    IMtr = IM[j - 2] + inctr

                    IM[j - 1] = IMtr
                    sf_x = IM[j - 1] / IM_geomean * self.g
                    sf_y = IM[j - 1] / IM_geomean * self.g
                    # Record into the outputs file
                    self.IM_output[rec, j-1] = IMtr
                    run = f"Record{rec + 1}_Run{j}"

                    # The trace intensity has been determined, now we can analyse
                    m = self.call_model()
                    self.time_series(dt, eq_name_x, eq_name_y, sf_x, sf_y)
                    if self.pflag:
                        print(f"[IDA] Record: {rec + 1}; Run: {j}; IM: {IMtr}")

                    th = SolutionAlgorithm(self.dt, dur, self.dcap, m.g.tnode, m.g.bnode, self.pflag, self.flag3d)
                    self.outputs[rec][j] = th.ntha_results

                    if th.c_index > 0:
                        # Stop tracing
                        tflag = 0
                        # Start filling
                        fflag = 1
                        # The value of j we traced to
                        jtrace = j
                        # Get the list of IMs
                        IMlist = IM
                        if j == jhunt:
                            print(f"[IDA] Warning: {run} - First trace for collapse resulted in collapse... ")
                    j += 1
                    op.wipe()

                # When the required resolution is reached, start filling until max_runs is reached
                # Make sure that number of runs are not exceeded
                if fflag == 1 and j <= self.max_runs:
                    print("[STEP] Filling the gaps...")

                    # Reorder the list so we can account for filled runs
                    IMlist = np.sort(IMlist)

                    # Determine the biggest gap in IM for the hunted runs
                    gap = 0.0
                    IMfil = 0.0
                    '''We go to the end of the list minus 1 because, if not we would be filling between a noncollapsing
                    and a collapsing run, for which we are not sure if that filling run would be a non collapse -
                    In short, does away with collapsing fills'''
                    for ii in range(1, len(IMlist) - 1):
                        # Find the running gap of hunted runs
                        temp = IMlist[ii] - IMlist[ii - 1]
                        if temp > gap:
                            # Update to maximum gap
                            gap = temp
                            # Determine new filling IM
                            IMfil = IMlist[ii - 1] + gap / 2

                    IM[j-1] = IMfil
                    IMlist = np.append(IMlist, IMfil)
                    sf_x = IM[j - 1] / IM_geomean * self.g
                    sf_y = IM[j - 1] / IM_geomean * self.g

                    # Record into the outputs file
                    self.IM_output[rec, j-1] = IMfil
                    run = f"Record{rec + 1}_Run{j}"

                    m = self.call_model()
                    self.time_series(dt, eq_name_x, eq_name_y, sf_x, sf_y)
                    if self.pflag:
                        print(f"[IDA] Record: {rec + 1}; Run: {j}; IM: {IMfil}")

                    th = SolutionAlgorithm(self.dt, dur, self.dcap, m.g.tnode, m.g.bnode, self.pflag, self.flag3d)
                    self.outputs[rec][j] = th.ntha_results

                    # Increment run number
                    j += 1
                    op.wipe()

                # Wrap it up and finish
                if j == self.max_runs and hflag == 1:
                    print('[IDA] Warning: Collapse not achieved, increase increment or number of runs...')
                if j == self.max_runs and fflag == 0:
                    print('[IDA] Warning: No filling, algorithm still tracing for collapse (reduce increment & '
                          'increase runs)...')
                op.wipe()

        print('[IDA] Finished IDA HTF')
