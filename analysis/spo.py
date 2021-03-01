"""
Performs static pushover analysis
"""
import openseespy.opensees as op
import numpy as np


class SPO:
    def __init__(self, cntr_node, disp_dir, base_cols, dref=1.0, nstep=1000, flag3d=False, direction=0):
        """
        Initialize static pushover definition
        :param cntr_node: int                       Node to control with displacement integrator
        :param disp_dir: int                        DOF the loading is applied
        :param base_cols: list(str)                 Base column IDs
        :param dref: float                          Reference displacement to which cycles are run. Corresponds to yield
                                                    or equivalent other, such as 1mm (in m)
        :param nstep: int                           Number of steps
        :param flag3d: bool                         True for 3D modelling, False for 2D modelling
        :param direction: int
        """
        self.TOL = 1e-08
        self.ITERINIT = 10
        self.dref = dref
        self.cntr_node = cntr_node
        self.disp_dir = disp_dir
        self.nstep = nstep
        self.base_cols = base_cols
        self.flag3d = flag3d
        self.direction = direction
        self.TEST_TYPE = 'NormDispIncr' if not self.flag3d else 'EnergyIncr'
        self.ALGORITHM_TYPE = 'KrylovNewton'
        self.NEGLIGIBLE = 1e-09

    def load_pattern(self, nodes, load_pattern=2, heights=None, mode_shape=None, nbays_x=None, nbays_y=None):
        """
        Define the load pattern
        :param nodes: list(int)                     Nodes (rightmost nodes) to which load pattern is applied to
        :param load_pattern: str                    Load pattern shape for static pushover analysis
                                                    0 = Uniform pattern
                                                    1 = Triangular pattern
                                                    2 = First-mode proportional pattern
        :param heights: list                        Storey heights (compatible with 1st load pattern)
        :param mode_shape: list                     1st mode shape (compatible with 2nd load pattern)
        :param nbays_x: int
        :param nbays_y: int
        :return: None
        """
        loads = []
        if load_pattern == 0:
            print('[STEP] Applying Uniform load pattern...')
            for i in nodes:
                loads.append(1.0)

        elif load_pattern == 1:
            print('[STEP] Applying triangular load pattern...')
            for h in range(len(heights)):
                if heights[h] != 0.0:
                    loads.append(heights[h] / heights[-1])

        elif load_pattern == 2:
            print('[STEP] Applying 1st mode proportional load pattern...')
            for i in mode_shape:
                loads.append(i)

        else:
            raise ValueError('[EXCEPTION] Wrong load pattern is supplied.')

        op.timeSeries('Linear', 4)
        op.pattern('Plain', 400, 4)
        if self.flag3d:
            # Number of stories
            nst = len(heights) - 1
            # Number of nodes
            n_nodes = (nbays_y + 1) * (nbays_x + 1)

            # Pushing all nodes with masses assigned to them
            for xbay in range(1, nbays_x + 2):
                for ybay in range(1, nbays_y + 2):
                    for st in range(1, nst + 1):
                        nodepush = int(f"{xbay}{ybay}{st}")
                        fpush = loads[st - 1]
                        if self.direction == 0:
                            op.load(nodepush, fpush, self.NEGLIGIBLE, self.NEGLIGIBLE, self.NEGLIGIBLE,
                                    self.NEGLIGIBLE, self.NEGLIGIBLE)
                        else:
                            op.load(nodepush, self.NEGLIGIBLE, fpush, self.NEGLIGIBLE, self.NEGLIGIBLE,
                                    self.NEGLIGIBLE, self.NEGLIGIBLE)

        else:
            for fpush, nodepush in zip(loads, nodes):
                if self.flag3d:
                    op.load(nodepush, fpush, self.NEGLIGIBLE, self.NEGLIGIBLE, self.NEGLIGIBLE, self.NEGLIGIBLE,
                            self.NEGLIGIBLE)
                else:
                    op.load(nodepush, fpush, self.NEGLIGIBLE, self.NEGLIGIBLE)

    def set_analysis(self, heights):
        """
        Sets up the initial analysis parameters
        :param heights: list
        :return: None
        """
        print('[INITIALIZE] Static pushover analysis has commenced...')
        if self.flag3d:
            op.constraints("Penalty", 1e15, 1e15)
            op.system('UmfPack')
        else:
            op.constraints('Plain')
            op.system('BandGeneral')
        op.numberer('RCM')
        op.test(self.TEST_TYPE, self.TOL, self.ITERINIT)
        op.algorithm(self.ALGORITHM_TYPE)
        op.integrator('DisplacementControl', self.cntr_node, self.disp_dir, 0.1 * heights[-1] / self.nstep)
        op.analysis('Static')

    def seek_solution(self):
        """
        Searches for a solution by using different test conditions or algorithms
        :return: ndarrays                           Top displacement vs Base shear
        """
        # It happens so, that column shear ID matches the disp_dir ID, they are not the same thing
        col_shear_idx = self.disp_dir

        '''Seek for a solution using different test conditions or algorithms'''
        # Set the initial values to start the while loop
        # The feature of disabling the possibility of having a negative loading has been included.
        #   adapted from a similar script by Prof. Garbaggio
        ok = 0
        step = 1
        loadf = 1.0

        # Recording top displacement and base shear
        topDisp = np.array([op.nodeResponse(self.cntr_node, self.disp_dir, 1)])
        baseShear = np.array([0.0])
        for col in self.base_cols:
            baseShear[0] += op.eleForce(int(col), col_shear_idx)

        while step <= self.nstep and ok == 0 and loadf > 0:
            ok = op.analyze(1)
            loadf = op.getTime()
            if ok != 0:
                print('[STEP] Trying relaxed convergence...')
                op.test(self.TEST_TYPE, self.TOL * 0.01, int(self.ITERINIT * 50))
                ok = op.analyze(1)
                op.test(self.TEST_TYPE, self.TOL, self.ITERINIT)
            if ok != 0:
                print('[STEP] Trying Newton with initial then current...')
                op.test(self.TEST_TYPE, self.TOL * 0.01, int(self.ITERINIT * 50))
                op.algorithm('Newton', '-initialThenCurrent')
                ok = op.analyze(1)
                op.algorithm(self.ALGORITHM_TYPE)
                op.test(self.TEST_TYPE, self.TOL, self.ITERINIT)
            if ok != 0:
                print('[STEP] Trying ModifiedNewton with initial...')
                op.test(self.TEST_TYPE, self.TOL * 0.01, int(self.ITERINIT * 50))
                op.algorithm('ModifiedNewton', '-initial')
                ok = op.analyze(1)
                op.algorithm(self.ALGORITHM_TYPE)
                op.test(self.TEST_TYPE, self.TOL, self.ITERINIT)
            if ok != 0:
                print('[STEP] Trying KrylovNewton...')
                op.test(self.TEST_TYPE, self.TOL * 0.01, int(self.ITERINIT * 50))
                op.algorithm('KrylovNewton')
                ok = op.analyze(1)
                op.algorithm(self.ALGORITHM_TYPE)
                op.test(self.TEST_TYPE, self.TOL, self.ITERINIT)
            if ok != 0:
                print('[STEP] Perform a Hail Mary...')
                op.test('FixedNumIter', self.ITERINIT)
                ok = op.analyze(1)

            # Recording the displacements and base shear forces
            topDisp = np.append(topDisp, op.nodeResponse(self.cntr_node, self.disp_dir, 1))
            eleForceTemp = 0.0
            for col in self.base_cols:
                eleForceTemp += op.eleForce(int(col), col_shear_idx)
            baseShear = np.append(baseShear, eleForceTemp)

            loadf = op.getTime()
            step += 1

        # Reverse sign of base_shear
        if min(baseShear) < 0.0:
            baseShear = -baseShear

        # Failure conditions
        if ok != 0:
            print('[FAILURE] DispControl analysis failed')
        else:
            print('[SUCCESS] DispControl analysis successful')

        if loadf <= 0:
            print(f"[FAILURE] Stopped because of load factor below zero: {loadf}")

        return topDisp, baseShear
