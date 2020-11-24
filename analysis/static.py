"""
Performs static analysis
"""
import openseespy.opensees as ops


class Static:
    def __init__(self):
        self.NSTEP = 10.0
        self.TOL = 1e-08

    def static_analysis(self):
        """
        Starts static analysis
        :return: None
        """
        dgravity = 1.0 / self.NSTEP
        ops.constraints('Transformation')
        ops.numberer('RCM')
        ops.system('UmfPack')
        ops.test('EnergyIncr', self.TOL, 10)
        ops.integrator('LoadControl', dgravity)
        ops.algorithm('Newton')
        ops.analysis('Static')
        ops.analyze(self.NSTEP)
        ops.loadConst('-time', 0.0)
