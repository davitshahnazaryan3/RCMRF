"""
Performs static analysis
"""
import openseespy.opensees as op


class Static:
    NSTEP = 1
    TOL = 1.e-08

    def static_analysis(self, path=None, flag3d=False):
        """
        Starts static analysis
        :param path: Path
        :param flag3d: bool
        :return: None
        """
        # Load increment
        dgravity = 1.0 / self.NSTEP
        # Determine next time step for an analysis
        op.integrator('LoadControl', dgravity)
        # Renumber dofs to minimize band-width (optimization)
        op.numberer('RCM')
        # Handling of boundary conditions
        if flag3d:
            op.constraints('Penalty', 1.0e15, 1.0e15)
            # Determine if convergence has been achieved at the end of an iteration step
            op.test('EnergyIncr', self.TOL, 10)
            # How to store and solve the system of equations in the analysis (large model: try UmfPack)
            op.system('UmfPack')
        else:
            op.constraints('Plain')
            op.test('NormDispIncr', self.TOL, 6)
            op.system('BandGeneral')

        # Use Newton's solution algorithm: updates tangent stiffness at every iteration
        op.algorithm('Newton')
        # Define type of analysis (static or transient)
        op.analysis('Static')
        # Apply gravity
        op.analyze(self.NSTEP)
        # Maintain constant gravity loads and reset time to zero
        op.loadConst('-time', 0.0)

        # Write to a tcl file
        if path:
            file = open(path / "Models/static.tcl", "w+")
            file.write("# Static analysis parameters")
            file.write("\nconstraints Penalty 1.0e15 1.0e15;")
            file.write("\nnumberer RCM;")
            file.write("\nsystem UmfPack;")
            file.write(f"\ntest NormDispIncr {self.TOL} 6;")
            file.write("\nalgorithm Newton;")
            file.write(f"\nintegrator LoadControl {dgravity};")
            file.write("\nanalysis Static;")
            file.write(f"\nanalyze {self.NSTEP};")
            file.write("\nloadConst -time 0.0;")
            file.close()
