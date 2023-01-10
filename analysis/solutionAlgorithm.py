"""
Performs nonlinear time history analysis
Displacements are in m
Accelerations are in g
Drifts are in %
NLRHA = Non-linear response history analysis
"""
import os
from scipy.interpolate import interp1d
import openseespy.opensees as op
import numpy as np
import warnings

from utils.utils import read_text, createFolder


class SolutionAlgorithm:
    g = 9.81
    ITER = 50
    ALGORITHM_TYPE = 'KrylovNewton'
    c_index = 0

    def __init__(self, dt, tmax, dc, tnode, bnode, dt_gm=None, eq_x=None, eq_y=None, f_x=1., f_y=1.,
                 pflag=True, flag3d=False, use_recorder=True, recorder_cache=None):
        """
        Procedure to execute the NLRHA
        :param dt: float                            Analysis time step
        :param tmax: float                          Length of the record (including padding of 0's)
        :param dc: float                            Drift capacity for both storey and roof drift (%)
        :param tnode: list(int)                     List of top nodes
        :param bnode: list(int)                     List of bottom nodes
        :param dt_gm: float
        :param eq_x: list
        :param eq_y: list
        :param f_x: float
        :param f_y: float
        :param pflag: bool                          Whether print information on screen or not
        :param flag3d: bool                         True for 3D modelling, False for 2D modelling
        :param use_recorder: bool                   Uses openseespy recorder to output file instead of node recorders
        :param recorder_cache: str                  Acceleration cache filename, don't leave empty when running multiple
                                                    analysis or MSA to avoid file rewrites
        """
        self.dt = dt
        self.tmax = tmax
        self.dc = dc
        self.dt_gm = dt_gm
        self.eq_x = eq_x
        self.eq_y = eq_y
        self.f_x = f_x
        self.f_y = f_y
        self.use_recorder = use_recorder
        self.pflag = pflag
        self.flag3d = flag3d
        self.recorder_cache = recorder_cache

        if self.flag3d:
            self.tnode = np.array(tnode)
            self.bnode = np.array(bnode)
        else:
            self.tnode = np.reshape(np.array(tnode), (1, len(tnode)))
            self.bnode = np.reshape(np.array(bnode), (1, len(bnode)))

        if not self.flag3d:
            self.TEST_TYPE = 'NormDispIncr'
            self.TOL = 1e-04
        else:
            self.TEST_TYPE = 'NormDispIncr'
            self.TOL = 1e-04
            # self.TEST_TYPE = 'EnergyIncr'
            # self.TOL = 1e-06

        # Interpolation functions
        if eq_x is not None and not use_recorder:
            # todo, not yet implemented
            self.int_x, self.int_y = self._create_interpolation_functions_for_accelerations()

        # Set analysis parameters
        self._set_analysis()
        # Run analysis
        self.ntha_results = self._seek_solution()

    def _create_interpolation_functions_for_accelerations(self):
        # Generate 0s to the end of eq_x and eq_y
        count_zero = int((self.tmax - len(self.eq_x) * self.dt_gm) / self.dt_gm)
        eq_x = self.eq_x * self.f_x
        eq_y = self.eq_y * self.f_y

        # eq_x = np.append(self.eq_x * self.f_x, np.zeros(count_zero))
        # eq_y = np.append(self.eq_y * self.f_x, np.zeros(count_zero))

        # Create time history
        time_history = np.linspace(self.dt_gm, round(self.tmax - 10 + self.dt_gm, 5),
                                   int(round((self.tmax - 10) / self.dt_gm, 0)))

        function_x = interp1d(time_history, eq_x, bounds_error=False, fill_value=0)
        function_y = interp1d(time_history, eq_y, bounds_error=False, fill_value=0)

        # try:
        #     # Create time history
        #     time_history = np.arange(self.dt_gm, self.tmax - 10 + self.dt_gm, self.dt_gm)
        #     # Create interpolation function
        #     function_x = interp1d(time_history, eq_x, bounds_error=False, fill_value=0)
        #     function_y = interp1d(time_history, eq_y, bounds_error=False, fill_value=0)
        # except:
        #     # Create time history
        #     time_history = np.arange(self.dt_gm, self.tmax - 10, self.dt_gm)
        #     # Create interpolation function
        #     function_x = interp1d(time_history, eq_x, bounds_error=False, fill_value=0)
        #     function_y = interp1d(time_history, eq_y, bounds_error=False, fill_value=0)

        return function_x, function_y

    def _set_analysis(self):
        """
        Sets up initial analysis parameters
        :return: None
        """
        op.test(self.TEST_TYPE, self.TOL, self.ITER)
        op.algorithm(self.ALGORITHM_TYPE)
        op.integrator('Newmark', 0.5, 0.25)
        op.analysis('Transient')

    def _calculate_ground_acceleration(self):
        """
        Calculates ground acceleration
        :return:
        """
        if not self.dt_gm:
            # Outputting relative accelerations
            return 0, 0

        # Otherwise calculate ground acceleration to add to relative acceleration amounting to absolute acceleration
        current_time = op.getTime()

        ground_acceleration_x = self.int_x(current_time)
        ground_acceleration_y = self.int_y(current_time)

        return ground_acceleration_x, ground_acceleration_y

    def _call_algorithms(self, ok, control_time):
        """
        Calls algorithms
        :param ok: int
        :param control_time: float
        :return: None
        """
        dtt = self.dt
        if ok != 0:
            if self.pflag:
                print(f"[FAILURE] Failed at {control_time} of {self.tmax} seconds")

        if self.eq_x is not None or self.use_recorder:
            if ok != 0:
                if self.pflag:
                    print(f"[FAILURE] Failed at {control_time} - Reduced timestep by half...")
                dtt = 0.5 * self.dt
                ok = op.analyze(1, dtt)
            if ok != 0:
                if self.pflag:
                    print(f"[FAILURE] Failed at {control_time} - Reduced timestep by quarter...")
                dtt = 0.25 * self.dt
                ok = op.analyze(1, dtt)

        if ok != 0:
            if self.pflag:
                print(f"[FAILURE] Failed at {control_time} - Trying Broyden...")
            op.algorithm('Broyden', 8)
            dtt = self.dt
            ok = op.analyze(1, dtt)
            op.algorithm(self.ALGORITHM_TYPE)
        if ok != 0:
            if self.pflag:
                print(f"[FAILURE] Failed at {control_time} - Trying Newton with initial tangent...")
            op.algorithm('Newton', '-initial')
            dtt = self.dt
            ok = op.analyze(1, dtt)
            op.algorithm(self.ALGORITHM_TYPE)
        if ok != 0:
            if self.pflag:
                print(f"[FAILURE] Failed at {control_time} - Trying NewtonWithLineSearch...")
            op.algorithm('NewtonLineSearch', 0.8)
            dtt = self.dt
            ok = op.analyze(1, self.dt)
            op.algorithm(self.ALGORITHM_TYPE)
        if ok != 0:
            if self.pflag:
                print(f"[FAILURE] Failed at {control_time} - Trying Newton with initial tangent & relaxed "
                      f"convergence...")
            op.test('NormDispIncr', self.TOL * 0.1, self.ITER * 50)
            op.algorithm('Newton', '-initial')
            dtt = self.dt
            ok = op.analyze(1, dtt)
            op.test(self.TEST_TYPE, self.TOL, self.ITER)
            op.algorithm(self.ALGORITHM_TYPE)
        if ok != 0:
            if self.pflag:
                print(f"[FAILURE] Failed at {control_time} - Trying NewtonWithLineSearch & relaxed convergence...")
            op.test('NormDispIncr', self.TOL * 0.1, self.ITER * 50)
            op.algorithm('NewtonLineSearch', 0.8)
            dtt = self.dt
            ok = op.analyze(1, dtt)
            op.test(self.TEST_TYPE, self.TOL, self.ITER)
            op.algorithm(self.ALGORITHM_TYPE)
        # Next, halve the timestep with both algorithm and tolerance reduction
        if self.eq_x is not None or self.use_recorder:
            if ok != 0:
                if self.pflag:
                    print(f"[FAILURE] Failed at {control_time} - Trying Newton with initial tangent, reduced timestep &"
                          f" relaxed convergence...")
                op.test('NormDispIncr', self.TOL * 0.1, self.ITER * 50)
                op.algorithm('Newton', '-initial')
                dtt = 0.5 * self.dt
                ok = op.analyze(1, dtt)
                op.test(self.TEST_TYPE, self.TOL, self.ITER)
                op.algorithm(self.ALGORITHM_TYPE)
            if ok != 0:
                if self.pflag:
                    print(f"[FAILURE] Failed at {control_time} - Trying NewtonWithLineSearch, reduced timestep &"
                          f" relaxed convergence...")
                op.test('NormDispIncr', self.TOL * 0.1, self.ITER * 50)
                op.algorithm('NewtonLineSearch', 0.8)
                dtt = 0.5 * self.dt
                ok = op.analyze(1, dtt)
                op.test(self.TEST_TYPE, self.TOL, self.ITER)
                op.algorithm(self.ALGORITHM_TYPE)

        if ok != 0:
            if self.pflag:
                print(f"[FAILURE] Failed at {control_time} - exit analysis...")
            self.c_index = -1

        self.dt_analysis = dtt

    def verify_against_zerolength(self, nst, h):
        for i in range(1, nst + 1):

            # Find the coordinates of the nodes in Global Y (2 for 2D, 3 for 3D)
            if self.flag3d:
                top2 = op.nodeCoord(int(self.tnode[0, i - 1]), 3)
                bot2 = op.nodeCoord(int(self.bnode[0, i - 1]), 3)
            else:
                top2 = op.nodeCoord(int(self.tnode[0, i - 1]), 2)
                bot2 = op.nodeCoord(int(self.bnode[0, i - 1]), 2)
            dist = top2 - bot2

            h = np.append(h, dist)

            if dist == 0:
                warnings.warn('[WARNING] Zerolength found in drift check.')
        return h

    def _seek_solution(self):
        """
        Performs the analysis and tries to find a converging solution
        :return: ndarrays                           Accelerations, displacements and storey drifts
        """
        createFolder("cache")

        # Set up analysis parameters
        control_time = 0.0
        ok = 0
        mdrift_init = 0.0

        # Number of dimensions to run analysis for
        d = 2 if self.flag3d else 1

        # Create recorders for each direction
        if self.use_recorder:
            for j in range(d):
                if self.recorder_cache is None:
                    self.recorder_cache = f'accelerations_{j + 1}.txt'

                op.recorder('Node', '-file', f'cache/{self.recorder_cache}',
                            '-timeSeries', int(f"5{j + 1}"),
                            '-nodeRange', int(self.bnode[0, 0]), int(self.tnode[0, -1]),
                            '-dof', j + 1, 'accel')

        # Set up the storey drift and acceleration values
        h = np.array([])
        # Number of storeys
        nst = self.tnode.shape[1]

        # initialize maximum drifts and accelerations
        mdrift = np.zeros((d, nst))
        maccel = np.zeros((d, nst + 1))

        # Recorders for displacements, accelerations and drifts, initialization
        displacements = np.zeros((d, nst + 1, 1))
        accelerations = np.zeros((d, nst + 1, 1))
        drifts = np.zeros((d, nst, 1))
        h = self.verify_against_zerolength(nst, h)

        # Run the actual analysis now
        while self.c_index == 0 and control_time <= self.tmax and ok == 0:
            # Start analysis
            ok = op.analyze(1, self.dt)
            control_time = op.getTime()

            # If the analysis fails, try the following changes to achieve convergence
            # Analysis will be slower in here though...
            self._call_algorithms(ok, control_time)

            # Recorders
            tempAccel = np.zeros((d, nst + 1, 1))
            tempDisp = np.zeros((d, nst + 1, 1))
            tempDrift = np.zeros((d, nst, 1))

            # Recording EDPs at each storey level to return
            # For each direction
            if not self.use_recorder:
                for j in range(tempAccel.shape[0]):
                    # At each storey level
                    for i in range(nst + 1):
                        if i == nst:
                            # Index 0 indicates along X direction, and 1 indicates along Y direction
                            # Nodal accelerations in g
                            tempAccel[j, i, 0] = op.nodeAccel(int(self.tnode[0, i - 1]), j + 1) / self.g
                            # Nodal displacements in m
                            tempDisp[j, i, 0] = op.nodeDisp(int(self.tnode[0, i - 1]), j + 1)
                        else:
                            # Get the PGA values (nodeAccel returns relative, not absolute values, so it will be 0)
                            tempAccel[j, i, 0] = op.nodeAccel(int(self.bnode[0, i]), j + 1) / self.g
                            tempDisp[j, i, 0] = op.nodeDisp(int(self.bnode[0, i]), j + 1)
                        if i > 0:
                            # Storey height
                            cht = h[i - 1]
                            # Storey drifts in %
                            tempDrift[j, i - 1, 0] = 100.0 * abs(tempDisp[j, i, 0] - tempDisp[j, i - 1, 0]) / cht

                # Appending into the global numpy arrays to return
                accelerations = np.append(accelerations, tempAccel, axis=2)

            displacements = np.append(displacements, tempDisp, axis=2)
            drifts = np.append(drifts, tempDrift, axis=2)

            # Check storey drifts and accelerations
            for i in range(1, nst + 1):
                # Top node displacement
                tnode_disp_x = op.nodeDisp(int(self.tnode[0, i - 1]), 1)
                tnode_disp_y = op.nodeDisp(int(self.tnode[0, i - 1]), 2)
                # Bottom node displacement
                bnode_disp_x = op.nodeDisp(int(self.bnode[0, i - 1]), 1)
                bnode_disp_y = op.nodeDisp(int(self.bnode[0, i - 1]), 2)
                # Storey height
                cht = h[i - 1]
                # Storey drift in %
                cdrift_x = 100.0 * abs(tnode_disp_x - bnode_disp_x) / cht
                cdrift_y = 100.0 * abs(tnode_disp_y - bnode_disp_y) / cht
                # Record the peak storey drifts
                if cdrift_x >= mdrift[0, i - 1]:
                    mdrift[0, i - 1] = cdrift_x

                if self.flag3d:
                    if cdrift_y >= mdrift[1, i - 1]:
                        mdrift[1, i - 1] = cdrift_y
                else:
                    cdrift_y = cdrift_x

                if max(cdrift_x, cdrift_y) >= mdrift_init:
                    mdrift_init = max(cdrift_x, cdrift_y)

            # Check whether drift capacity has been exceeded
            if mdrift_init >= self.dc:
                # If it was exceeded then local structure collapse is assumed
                self.c_index = 1
                # Hard cap the mdrift_init value to the drift capacity
                mdrift_init = self.dc

        # Wipe the model
        op.wipe()

        # Record the absolute accelerations, when use_recorder is True
        if self.use_recorder:
            accelerations = []
            for j in range(d):
                if self.recorder_cache is None:
                    self.recorder_cache = f'accelerations_{j + 1}.txt'

                accelerations.append(np.transpose(read_text(f'cache/{self.recorder_cache}') / self.g))

                os.remove(f'cache/{self.recorder_cache}')

        accelerations = np.array(accelerations)

        if self.c_index == -1:
            print(f"[FAILURE] Analysis failed to converge at {control_time} of {self.tmax}.")
        if self.c_index == 0:
            print('[SUCCESS] Analysis completed successfully.')
        if self.c_index == 1:
            print('[FAILURE] Local structure collapse.')

        return accelerations, displacements, drifts
