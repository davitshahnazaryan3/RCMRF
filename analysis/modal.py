"""
Performs eigenvalue analysis
"""
import openseespy.opensees as op
import numpy as np


class Modal:
    def __init__(self, num_modes, damp_modes=None, damping=0.05, path=None):
        """
        Initializes modal analysis
        :param num_modes: int                   Number of modes of interest
        :param damp_modes: list(int)            2 element List of damping modes (e.g. [1, 3])
        :param damping: float                   Ratio of critical damping to be applied to the listed modes
        :param path: bool                       Exporting Model to path
        """
        self.num_modes = num_modes
        self.path = path
        self.damp_modes = damp_modes
        self.damping = damping
        self.lam = self._compute_eigenvectors()
        self._record_stuff()
        self.omega, self.freq, self.period = self._extract_eigenvalues(self.lam)

        self.xi_modes = self._get_damping(self.omega)

        self.file = None

    def _compute_eigenvectors(self):
        """
        Computes eigen values
        :return: float                          Eigenvalue
        """

        lam = None
        try:
            lam = op.eigen(self.num_modes)
        except:
            print('[EXCEPTION] Eigensolver failed, trying genBandArpack...')
            try:
                lam = op.eigen('-genBandArpack', self.num_modes)
            except:
                print('[EXCEPTION] Eigensolver failed, trying fullGenLapack...')
                try:
                    lam = op.eigen('-fullGenLapack', self.num_modes)
                except:
                    print('[EXCEPTION] Eigensolver failed, trying symmBandLapack...')
                    try:
                        lam = op.eigen('-symmBandLapack', self.num_modes)
                    except:
                        print('[EXCEPTION] Eigensolver failed.')

        # Write to file
        if self.path:
            self.file = open(self.path / "Models/modal_analysis.tcl", "w+")
            self.file.write("# Modal analysis procedure")
            self.file.write("\n# Solve for lambda")
            self.file.write(f"\nset lambda [eigen {self.num_modes}];")

            self.file.write("\n# If solver failed, try another")
            self.file.write("\nif {[lindex $lambda 0] <= 0} {")
            self.file.write(f"\n\tset lambda [eigen -genBandArpack  {self.num_modes}];")
            self.file.write('\n\tputs "Eigensolver failed, trying genBandArpack...";\n}')

            self.file.write("\nif {[lindex $lambda 0] <= 0} {")
            self.file.write(f"\n\tset lambda [eigen -fullGenLapack  {self.num_modes}];")
            self.file.write('\n\tputs "Eigensolver failed, trying fullGenLapack...";\n}')

            self.file.write("\nif {[lindex $lambda 0] <= 0} {")
            self.file.write(f"\n\tset lambda [eigen -symmBandLapack  {self.num_modes}];")
            self.file.write('\n\tputs "Eigensolver failed, trying symmBandLapack...";\n}')

            self.file.write("\nif {[lindex $lambda 0] <= 0} {")
            self.file.write('\n\tputs "Eigensolver failed.";\n}')

            self.file.write("\n\n# Record the eigenvectors")
            self.file.write("\nrecord")

        return lam

    @staticmethod
    def _record_stuff():
        """
        Records the eigenvectors
        :return: None
        """
        op.record()

    def _extract_eigenvalues(self, lam):
        """
        Extracts eigenvalues to appropriate arrays
        :param lam: float                       Eigenvalue
        :return: lists                          Circular frequencies, frequencies and periods
        """
        omega = []
        freq = []
        period = []
        for m in range(self.num_modes):
            omega.append(np.sqrt(lam[m]))
            freq.append(np.sqrt(lam[m]) / 2 / np.pi)
            period.append(2 * np.pi / np.sqrt(lam[m]))

        # write to file
        if self.path:
            self.file.write("\n\n# Extract the eigenvalues to the appropriate arrays")
            self.file.write("\nset omega {};")
            self.file.write("\nset freq {};")
            self.file.write("\nset periods {};")

            self.file.write("\nforeach lam $lambda {")
            self.file.write("\n\tlappend omega [expr sqrt($lam)]")
            self.file.write("\n\tlappend freq [expr sqrt($lam)/(2*3.14159)]")
            self.file.write("\n\tlappend periods [expr (2*3.14159)/sqrt($lam)]\n};")

            self.file.close()

        return omega, freq, period

    def _get_damping(self, omega):
        """
        Computes Rayleigh damping
        :param omega: list                      List of circular frequencies
        :return: list                           Rayleigh damping
        """
        if self.damp_modes is None:
            self.damp_modes = [0]

        n = len(self.damp_modes)
        if n > 1:
            wi = omega[(self.damp_modes[0] - 1)]
            wj = omega[(self.damp_modes[1] - 1)]
            a0 = 2 * wi * wj / (wj ** 2 - wi ** 2) * (wj * self.damping - wi * self.damping)
            a1 = 2 * wi * wj / (wj ** 2 - wi ** 2) * (-self.damping / wj + self.damping / wi)
            modes = []
            for m in range(self.num_modes):
                wn = omega[m]
                modes.append(0.5 * (a0 / wn + a1 * wn))
            else:
                return modes

        elif n == 1:
            modes = 0.0
            return modes

        else:
            raise TypeError('[EXCEPTION] No damping mode was provided')
