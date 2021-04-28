"""
Performs eigenvalue analysis
"""
import openseespy.opensees as op
import numpy as np


class Modal:
    def __init__(self, num_modes, damp_modes=None, damping=0.05):
        """
        Initializes modal analysis
        :param num_modes: int                   Number of modes of interest
        :param damp_modes: list(int)            2 element List of damping modes (e.g. [1, 3])
        :param damping: float                   Ratio of critical damping to be applied to the listed modes
        """
        self.num_modes = num_modes
        self.damp_modes = damp_modes
        self.damping = damping
        self.lam = self.compute_eigenvectors()
        self.record_stuff()
        self.omega, self.freq, self.period = self.extract_eigenvalues(self.lam)

        self.xi_modes = self.get_damping(self.omega)

    def compute_eigenvectors(self):
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

        return lam

    @staticmethod
    def record_stuff():
        """
        Records the eigenvectors
        :return:
        """
        op.record()

    def extract_eigenvalues(self, lam):
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

        return omega, freq, period

    def get_damping(self, omega):
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
