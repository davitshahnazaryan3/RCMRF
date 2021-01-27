"""
Performs nonlinear time history analysis
Displacements are in m
Accelerations are in g
Drifts are in %
"""
import openseespy.opensees as op
import numpy as np, warnings


class SolutionAlgorithm:
    def __init__(self, dt, tmax, dc, tnode, bnode, pflag=True, flag3d=False):
        """
        Procedure to execute the NRHA of a 2D model
        :param dt: float                            Analysis time step
        :param tmax: float                          Length of the record (including padding of 0's)
        :param dc: float                            Drift capacity for both storey and roof drift (%)
        :param tnode: list(int)                     List of top nodes
        :param bnode: list(int)                     List of bottom nodes
        :param pflag: bool                          Whether print information on screen or not
        :param flag3d: bool                         True for 3D modelling, False for 2D modelling
        """
        self.dt = dt
        self.tmax = tmax
        self.dc = dc
        self.tnode = np.array(tnode)
        self.bnode = np.array(bnode)
        self.pflag = pflag
        self.flag3d = flag3d
        self.TEST_TYPE = 'NormDispIncr'
        self.TOL = 1e-04
        self.ITER = 50
        self.ALGORITHM_TYPE = 'KrylovNewton'
        self.c_index = 0
        self.set_analysis()
        self.ntha_results = self.seek_solution()

    def set_analysis(self):
        """
        Sets up initial analysis parameters
        :return: None
        """
        op.test(self.TEST_TYPE, self.TOL, self.ITER)
        op.algorithm(self.ALGORITHM_TYPE)
        op.integrator('Newmark', 0.5, 0.25)
        op.analysis('Transient')

    def seek_solution(self):
        """
        Performs the analysis and tries to find a converging solution
        :return: ndarrays                           Accelerations, displacements and storey drifts
        """
        # Set up analysis parameters
        control_time = 0.0
        ok = 0
        mdrift_init = 0.0

        # Set up the storey drift and acceleration values
        h = np.array([])
        nst = len(self.tnode)
        # maccel = np.zeros((nst+1, ))
        mdrift = np.zeros((nst, ))

        # Recorders for displacements, accelerations and drifts
        displacements = np.zeros((nst + 1, 1))
        accelerations = np.zeros((nst + 1, 1))
        drifts = np.zeros((nst, 1))
        bdg_h = 0.0
        for i in range(1, nst + 1):
            # Find the coordinates of the nodes in Global Y (2 for 2D, 3 for 3D)
            if self.flag3d:
                top2 = op.nodeCoord(int(self.tnode[i - 1]), 3)
                bot2 = op.nodeCoord(int(self.bnode[i - 1]), 3)
            else:
                top2 = op.nodeCoord(int(self.tnode[i - 1]), 2)
                bot2 = op.nodeCoord(int(self.bnode[i - 1]), 2)
            dist = top2 - bot2

            bdg_h = bdg_h + dist
            h = np.append(h, dist)

            if dist == 0:
                warnings.warn('[WARNING] Zerolength found in drift check.')

        # Run the actual analysis now
        while self.c_index == 0 and control_time <= self.tmax and ok == 0:
            # Start analysis
            ok = op.analyze(1, self.dt)
            control_time = op.getTime()

            # If the analysis fails, try the following changes to achieve convergence
            # Analysis will be slower in here though...

            if ok != 0:
                if self.pflag:
                    print(f"[FAILURE] Failed at {control_time} of {self.tmax} seconds")

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
                ok = op.analyze(1, self.dt)
                op.algorithm(self.ALGORITHM_TYPE)
            if ok != 0:
                if self.pflag:
                    print(f"[FAILURE] Failed at {control_time} - Trying Newton with initial tangent...")
                op.algorithm('Newton', '-initial')
                ok = op.analyze(1, self.dt)
                op.algorithm(self.ALGORITHM_TYPE)
            if ok != 0:
                if self.pflag:
                    print(f"[FAILURE] Failed at {control_time} - Trying NewtonWithLineSearch...")
                op.algorithm('NewtonLineSearch', 0.8)
                ok = op.analyze(1, self.dt)
                op.algorithm(self.ALGORITHM_TYPE)
            if ok != 0:
                if self.pflag:
                    print(f"[FAILURE] Failed at {control_time} - Trying Newton with initial tangent & relaxed "
                          f"convergence...")
                op.test('NormDispIncr', self.TOL * 0.1, self.ITER * 50)
                op.algorithm('Newton', '-initial')
                ok = op.analyze(1, self.dt)
                op.test(self.TEST_TYPE, self.TOL, self.ITER)
                op.algorithm(self.ALGORITHM_TYPE)
            if ok != 0:
                if self.pflag:
                    print(f"[FAILURE] Failed at {control_time} - Trying NewtonWithLineSearch & relaxed convergence...")
                op.test('NormDispIncr', self.TOL * 0.1, self.ITER * 50)
                op.algorithm('NewtonLineSearch', 0.8)
                ok = op.analyze(1, self.dt)
                op.test(self.TEST_TYPE, self.TOL, self.ITER)
                op.algorithm(self.ALGORITHM_TYPE)

            # Next, halve the timestep with both algorithm and tolerance reduction
            if ok != 0:
                if self.pflag:
                    print(f"[FAILURE] Failed at {control_time} - Trying Newton with initial tangent, reduced timestep &"
                          f" relaxed convergence...")
                op.test('NormDispIncr', self.TOL * 0.1, self.ITER * 50)
                op.algorithm('Newton', '-initial')
                ok = op.analyze(1, 0.5 * self.dt)
                op.test(self.TEST_TYPE, self.TOL, self.ITER)
                op.algorithm(self.ALGORITHM_TYPE)
            if ok != 0:
                if self.pflag:
                    print(f"[FAILURE] Failed at {control_time} - Trying NewtonWithLineSearch, reduced timestep &"
                          f" relaxed convergence...")
                op.test('NormDispIncr', self.TOL * 0.1, self.ITER * 50)
                op.algorithm('NewtonLineSearch', 0.8)
                ok = op.analyze(1, 0.5 * self.dt)
                op.test(self.TEST_TYPE, self.TOL, self.ITER)
                op.algorithm(self.ALGORITHM_TYPE)
            if ok != 0:
                if self.pflag:
                    print(f"[FAILURE] Failed at {control_time} - exit analysis...")
                self.c_index = -1

            # Recorders
            # TODO, modify for 3D flag, add recorders for Y direction, as the third axis in the following recorders
            tempAccel = np.zeros((nst + 1, 1))
            tempDisp = np.zeros((nst + 1, 1))
            tempDrift = np.zeros((nst, 1))

            # Recording EDPs at each storey level to return
            for i in range(nst + 1):
                if i == nst:
                    # TODO, add option for recorders to create separate txt files and a recorder with timeseries
                    #  for accelerations
                    # Nodal accelerations in g
                    tempAccel[i] = op.nodeAccel(int(self.tnode[i - 1]), 1) / 9.81
                    # Nodal displacements in m
                    tempDisp[i] = op.nodeDisp(int(self.tnode[i - 1]), 1)
                else:
                    tempAccel[i] = op.nodeAccel(int(self.bnode[i]), 1) / 9.81
                    tempDisp[i] = op.nodeDisp(int(self.bnode[i]), 1)
                if i > 0:
                    # Storey height
                    cht = h[i - 1]
                    # Storey drifts in %
                    tempDrift[i - 1] = 100.0 * abs(tempDisp[i] - tempDisp[i - 1]) / cht

            # Appending into the global numpy arrays to return
            accelerations = np.append(accelerations, tempAccel, axis=1)
            displacements = np.append(displacements, tempDisp, axis=1)
            drifts = np.append(drifts, tempDrift, axis=1)

            # # Base accelerations in [g] (this is going to be zero...)
            # base_accel = op.nodeAccel(int(self.bnode[0]), 1) / 9.81
            # if base_accel > maccel[0]:
            #     maccel[0] = base_accel

            # Check storey drifts and accelerations
            for i in range(1, nst+1):
                # Top node displacement
                tnode_disp = op.nodeDisp(int(self.tnode[i - 1]), 1)
                # Bottom node displacement
                bnode_disp = op.nodeDisp(int(self.bnode[i - 1]), 1)
                # Storey height
                cht = h[i - 1]
                # Storey drift in %
                cdrift = 100.0 * abs(tnode_disp - bnode_disp) / cht
                # Record the peak storey drifts
                if cdrift >= mdrift[i - 1]:
                    mdrift[i - 1] = cdrift

                if cdrift >= mdrift_init:
                    mdrift_init = cdrift

                # # Accelerations in g (this is not quite important)
                # node_accel = op.nodeAccel(int(self.tnode[i - 1]), 1) / 9.81
                # caccel = 1.0 * abs(node_accel)
                # if caccel >= maccel[i]:
                #     maccel[i] = caccel

            # Check whether drift capacity has been exceeded
            if mdrift_init >= self.dc:
                # If it was exceeded then local structure collapse is assumed
                self.c_index = 1
                # Hard cap the mdrift_init value to the drift capacity
                mdrift_init = self.dc
                # Wipe the model
                op.wipe()

        if self.c_index == -1:
            print(f"[FAILURE] Analysis failed to converge at {control_time} of {self.tmax}.")
        if self.c_index == 0:
            print('[SUCCESS] Analysis completed successfully.')
        if self.c_index == 1:
            print('[FAILURE] Local structure collapse.')

        return accelerations, displacements, drifts
