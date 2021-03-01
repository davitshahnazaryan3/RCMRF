"""
Performs static analysis
"""
import openseespy.opensees as op


class Static:
    def __init__(self):
        self.NSTEP = 10
        self.TOL = 1e-08

    def static_analysis(self, flag3d=False):
        """
        Starts static analysis
        :param flag3d: bool
        :return: None
        """
        # Load increment
        dgravity = 1.0 / self.NSTEP
        # Determine next time step for an analysis
        op.integrator('LoadControl', dgravity)
        # Renumber dofs to minimize band-width (optimization)
        op.numberer('RCM')
        # How to store and solve the system of equations in the analysis (large model: try UmfPack)
        op.system('UmfPack')
        # Handling of boundary conditions
        if flag3d:
            op.constraints('Penalty', 1.0e15, 1.0e15)
        else:
            op.constraints('Transformation')
        # Determine if convergence has been achieved at the end of an iteration step
        op.test('EnergyIncr', self.TOL, 10)
        # Use Newton's solution algorithm: updates tangent stiffness at every iteration
        op.algorithm('Newton')
        # Define type of analysis (static or transient)
        op.analysis('Static')
        # Apply gravity
        op.analyze(self.NSTEP)
        # Maintain constant gravity loads and reset time to zero
        op.loadConst('-time', 0.0)

