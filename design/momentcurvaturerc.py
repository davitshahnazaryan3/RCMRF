"""
A software to obtain Moment-curvature relationship of an RC element
Follows the recommendations of Prestressed concrete structures by Collins
Written to optimize for a reinforcement by Davit Shahnazaryan

Units
m for length
kN for forces
mpa for stress
m2 for area
1/m for curvature
+ for tension
- for compression
"""
import numpy as np
from scipy import optimize
import math

from design.plasticity import Plasticity
from utils.utils import getIndex
import warnings

warnings.filterwarnings('ignore')


class MomentCurvatureRC:
    def __init__(self, b, h, m_target, length=0., nlayers=0, p=0., d=.03, fc_prime=25, fy=415, young_mod_s=200e3,
                 soft_method="Collins", k_hard=1.0, fstiff=0.5, AsTotal=None, distAs=None):
        """
        init Moment curvature tool
        :param b: float                         Element sectional width
        :param h: float                         Element sectional height
        :param m_target: float                  Target flexural capacity
        :param length: float                    Distance from critical section to point of contraflexure
        :param nlayers: int                     Number of flexural reinforcement layers
        :param p: float                         Axial load (negative=compression, positive=tension)
        :param d: float                         Flexural reinforcement cover in m
        :param fc_prime: float                  Concrete compressive strength
        :param fy: float                        Reinforcement yield strength
        :param young_mod_s: float               Young modulus of reinforcement
        :param soft_method: str                 Method for the softening slope calculation
        :param k_hard: float                    Hardening slope of reinforcement (i.e. fu/fy)
        :param fstiff: float                    Stiffness reduction factor (50% per Eurocode 8), for the model only
        :param AsTotal: float                   Total reinforcement area (for beams only)
        :param distAs: list                     Relative distributions of reinforcement (for beams only)
        """
        self.b = b
        self.h = h
        self.m_target = m_target
        self.length = length
        self.nlayers = nlayers
        self.p = p
        self.d = d
        self.fc_prime = fc_prime
        self.fy = fy
        self.young_mod_s = young_mod_s
        self.EPSSH = 0.008
        self.EPSUK = 0.075
        self.k_hard = k_hard
        self.soft_method = soft_method
        self.fstiff = fstiff
        self.mi = np.nan
        self.epss = np.nan
        self.fst = np.nan
        self.phii = np.nan
        self.AsTotal = AsTotal
        self.distAs = distAs

        if self.distAs is None:
            self.distAs = np.array([0.5, 0.5])
        # Transverse reinforcement spacing in [m]
        self.TRANSVERSE_SPACING = 0.1
        # Transverse reinforcement diameter in [m]
        self.TRANSVERSE_DIAMETER = 8 / 1000
        # Number of transverse reinforcement legs
        self.TRANSVERSE_LEGS = 4
        # Residual strength ratio
        self.RESIDUAL = 0.0

    def _get_reinforcement_information(self, rebar):
        if self.nlayers == 0:
            z = np.array([self.h - self.d, self.d])
            rebar = np.array([rebar * self.distAs[1], rebar * self.distAs[0]])
        elif self.nlayers == 1:
            z = np.array([self.h - self.d, self.h / 2, self.d])
            rebar = np.array([rebar * 3 / 8, rebar * 2 / 8, rebar * 3 / 8])
        elif self.nlayers == 2:
            z = np.array([self.h - self.d, ((self.h - 2 * self.d) / (1 + self.nlayers) + self.d) * 2 - self.d,
                          ((self.h - 2 * self.d) / (1 + self.nlayers) + self.d), self.d])
            rebar = np.array([rebar * 4 / 12, rebar * 2 / 12, rebar * 2 / 12, rebar * 4 / 12])
        else:
            raise ValueError(f"[EXCEPTION] Wrong input: number of reinforcement layers is {self.nlayers}, should be 2 "
                             f"or lower!")

        return z, rebar

    def compute_stress(self, epsc, epsc_prime, c, z, residual=False):
        # reinforcement properties
        ey = self.fy / self.young_mod_s
        fu = self.k_hard * self.fy

        # Block parameters
        b1 = (4 - epsc / epsc_prime) / (6 - 2 * epsc / epsc_prime)
        a1b1 = epsc / epsc_prime - 1 / 3 * (epsc / epsc_prime) ** 2
        if epsc > 2 * epsc_prime:
            b1 = (4 - 2 * epsc_prime / epsc_prime) / (6 - 4 * epsc_prime / epsc_prime)
            a1b1 = 2 * epsc_prime / epsc_prime - 1 / 3 * (2 * epsc_prime / epsc_prime) ** 2

        # Top rebar strain
        if not residual:
            epss = (c - (self.h - z)) / c * epsc
        else:
            epss = (c - (self.h - z[:-1])) / c * epsc

        # Stresses
        stress = np.zeros(len(epss))
        for i in range(len(stress)):
            if abs(epss[i]) <= ey:
                stress[i] = self.young_mod_s * epss[i]
            elif ey < abs(epss[i]) <= self.EPSSH:
                stress[i] = np.sign(epss[i]) * self.fy
            elif self.EPSSH < abs(epss[i]) <= self.EPSUK:
                stress[i] = np.sign(epss[i]) * self.fy + (np.sign(epss[i]) * fu - np.sign(epss[i]) * self.fy) * \
                            np.sqrt((epss[i] - np.sign(epss[i]) * self.EPSSH) / (
                                    np.sign(epss[i]) * self.EPSUK - np.sign(epss[i]) * self.EPSSH))
            else:
                # This is an approximation, generally reinforcement will be sufficient enough not to surpass ultimate
                # strain. However, for calculation purpose this will be left here for now
                stress[i] = np.sign(epss[i]) * self.fy + (np.sign(epss[i]) * fu - np.sign(epss[i]) * self.fy) * \
                            np.sqrt((epss[i] - 2 * (np.sign(epss[i]) * self.EPSUK - np.sign(epss[i]) * self.EPSSH)) / (
                                    np.sign(epss[i]) * self.EPSUK - 2 * (np.sign(epss[i]) * self.EPSUK
                                                                         - np.sign(epss[i]) * self.EPSSH)))

        # Internal force in compressed concrete
        cc = c * a1b1 * self.fc_prime * self.b * 1000

        # Compressed height
        compr_height = min(b1 * c, self.h)
        return cc, compr_height, stress, epss

    def objective(self, c, data):
        """
        Objective function solving for compressed concrete height
        :param c: numpy.ndarray                 Compressed concrete height
        :param data: list                       Reinforcement characteristics
        :return: float                          Difference between internal and analysis forces
        """
        # Force it to look for only positive values of c
        c = abs(c)
        # Concrete strains
        epsc = data[0]
        epsc_prime = data[1]
        # Reinforcement
        rebar = data[2]

        # z = locations of rebar layers, from top to bottom (bottom corner = 0.0)
        z, rebar = self._get_reinforcement_information(rebar)

        # Compute stresses and internal forces
        cc, compr_height, stress, epss = self.compute_stress(epsc, epsc_prime, c, z)

        # Forces
        nslist = rebar * stress * 1000
        nint = cc + sum(nslist)

        self.mi = cc * (self.h / 2 - compr_height / 2) + sum(nslist * (z - self.h / 2))
        self.fst = abs(stress[-1])
        self.epss = abs(epss[-1])
        self.phii = epsc / c

        return abs(nint + self.p)

    def get_residual_strength(self, c, data):
        """
        Objective function solving for compressed concrete height
        Only top reinforcement is working under tension, bottom reinforcement is assumed to be lost
        Hence, compressed concrete height is within the cover of the top reinforcement
        :param c: numpy.ndarray                 Compressed concrete height
        :param data: list                       Reinforcement characteristics
        :return: float                          Difference between internal and analysis forces
        """
        # Force it to look for only positive values of c
        c = abs(c)
        # Concrete strains
        epsBot = data[0]
        epsc_prime = data[1]
        # Reinforcement
        rebar = data[2]

        # z = locations of rebar layers, from top to bottom (bottom corner = 0.0)
        z, rebar = self._get_reinforcement_information(rebar)

        # Get strain at top concrete fiber
        epsc = float(c * epsBot / (z[0] - c))

        # Compute stresses and internal forces
        cc, compr_height, stress, _ = self.compute_stress(epsc, epsc_prime, c, z, residual=True)

        # Internal forces in reinforcement bars
        ns = rebar[:-1] * stress * 1000

        # Total internal force, concrete + reinforcement
        nint = cc + sum(ns)

        # Internal bending moment and curvature
        self.mi = float(cc * (self.h / 2 - compr_height / 2) + sum(ns * (z[:-1] - self.h / 2)))
        self.phii = float(epsc / c)

        return abs(nint + self.p)

    def find_fracturing(self, epsc_prime, asinit):
        """
        Find fracturing point
        :param epsc_prime: float                    Concrete strain at max compressive strength
        :param asinit: float                        Total reinforcement area
        :return: None
        """
        # Take 10% more than the last recorded bottom reinforcement strain
        epss_bot = 1.1 * self.epss
        # epss_bot = 0.044
        # Initialize moment and compressed concrete height
        c = 0.01
        c = float(optimize.fsolve(self.get_residual_strength, c, [epss_bot, epsc_prime, asinit], factor=0.1))

        moment = self.mi
        phii = self.phii

        if moment >= 0:
            record = [phii, moment]
        else:
            record = [phii, 0.0]

        return record

    def max_moment(self, asi, epsc_prime):
        """
        Gets the maximum moment capacity
        :param asi: numpy.ndarray                   Total reinforcement area
        :param epsc_prime: float                    Concrete strain at peak compressive strength
        :return: float                              Maximum moment of element
        """
        asinit = asi[0]
        c = np.array([0.05])
        c = float(abs(optimize.fsolve(self.objective, c, [2 * epsc_prime, epsc_prime, asinit], factor=0.1)))
        return abs(self.mi / self.k_hard - self.m_target)

    def get_softening_slope(self, **kwargs):
        """
        defines the softening slope of the moment-curvature relationship
        :param kwargs: floats                       Total reinforcement area, curvature at yield, axial load ratio
                                                    transverse steel ratio
        :return: float                              Curvature and Moment at 0 and plastic hinge length
        """
        lp = Plasticity(lp_name="Priestley", db=20, fy=self.fy, fu=self.fy * self.k_hard, lc=self.length).get_lp()
        if self.soft_method == "Haselton":
            phiy = kwargs.get('curvature_yield', None)
            mu_phi = kwargs.get('curvature_ductility', None)
            nu = kwargs.get('axial_load_ratio', 0)
            ro_sh = kwargs.get('transverse_steel_ratio', None)
            if ro_sh is not None:
                theta_pc = min(0.76 * 0.031 ** nu * (0.02 + 40 * ro_sh) ** 1.02, 0.1)
            else:
                theta_pc = 0.1
            # TODO, fix phi_pc formula, the accuracy needs to be increased as it does not account for elastic portion
            phi_pc = theta_pc / lp
            phi_critical = phiy * mu_phi + phi_pc
            m_critical = 0.

        elif self.soft_method == "Collins":

            rebar_area = kwargs.get('rebar_area', None)
            young_modulus_rc = (3320 * np.sqrt(self.fc_prime) + 6900)
            n = .8 + self.fc_prime / 17
            epsc_prime = self.fc_prime / young_modulus_rc * n / (n - 1)
            record = self.find_fracturing(epsc_prime, rebar_area)

            # Fracturing curvature
            m_critical = record[1]
            phi_critical = record[0]

        else:
            raise ValueError("[EXCEPTION] Wrong method for the definition of softening slope!")
        return phi_critical, m_critical, lp

    def get_mphi(self, check_reinforcement=False, reinf_test=0., m_target=None, reinforcements=None, cover=None):
        # TODO, a bit too rigid, make it more flexible, easier to manipulate within IPBSD to achieve optimized designs
        # TODO, issue where fracturing curvature is not computed correctly and is equal to hardening curvature,
        #  look into it
        """
        Gives the Moment-curvature relationship
        :param check_reinforcement: bool            Gets moment for reinforcement provided (True) or applied
                                                    optimization for Mtarget (False)
        :param reinf_test: int                      Reinforcement for test
        :param m_target: float                      Target bending moment. This is a value that may be increased
                                                    depending on local ductility requirements
        :param reinforcements: list                 Positive and negative reinforcements (for beams only)
        :param cover: float                         Reinforcement cover, generally input when warnMin was triggered
        :return: dict                               M-phi response data, reinforcement and concrete data for detailing
        """
        if reinforcements is not None:
            reinforcements = np.array(reinforcements)
            self.AsTotal = sum(reinforcements)
            self.distAs = reinforcements / self.AsTotal
        if m_target is not None:
            self.m_target = float(m_target)
        if cover is not None:
            self.d = cover

        # Concrete properties
        # Assumption - parabolic stress-strain relationship for the concrete
        # concrete elasticity modulus MPa
        young_modulus_rc = (3320 * np.sqrt(self.fc_prime) + 6900)
        # young_modulus_rc = 22*((fc_prime+8)/10)**.3*1000
        n = .8 + self.fc_prime / 17
        k_parameter = 0.67 + self.fc_prime / 62
        epsc_prime = self.fc_prime / young_modulus_rc * n / (n - 1)
        # Reinforcement properties (500C grade)
        ey = self.fy / self.young_mod_s
        area = self.h * self.b
        inertia = self.b * self.h ** 3 / 12
        # Cracking moment calculation (irrelevant for the design, but will store the data for possible checks)
        lam_nw = 1  # for normal weight concrete
        fcr = 0.33 * lam_nw * np.sqrt(self.fc_prime)
        m_cr = (-self.p / area + fcr * 1000) * inertia / (self.h / 2)
        epscr = fcr / young_modulus_rc
        fcr_t = (-self.p / area + m_cr * self.h / 2 / inertia) / 1000
        yc = fcr * self.h / (fcr + fcr_t)
        phi_cr = epscr / yc

        ''' The "Process" '''
        epsc = np.linspace(epsc_prime * 2 / 500, 10 * epsc_prime, 400)
        sigma_c = self.fc_prime * n * epsc / epsc_prime / (n - 1 + np.power(epsc / epsc_prime, n * k_parameter))
        m = np.array([0])
        phi = np.array([0])
        sigmat = np.array([0])
        eps_tensile = np.array([0])

        # Optimize for longitudinal reinforcement at peak capacity
        if self.AsTotal is not None:
            asinit = np.array([self.AsTotal])
        else:
            asinit = np.array([0.005])

        asinit = abs(float(optimize.fsolve(self.max_moment, asinit, epsc_prime, factor=0.1)))

        # Are we doing a reinforcement check? If, yes...
        if check_reinforcement:
            c = np.array([0.05])
            self.mi = None
            init_factor = 2.
            while self.mi is None or np.isnan(self.mi[0]):
                c = abs(float(optimize.fsolve(self.objective, c, [init_factor * epsc_prime, epsc_prime, reinf_test],
                                              factor=0.1)))
                init_factor -= 0.1
            return self.mi

        # If not, get the full M-Phi curve
        else:
            for i in range(len(epsc)):
                # compressed section height optimization - make a good guess, otherwise convergence won't be achieved
                c = 0.05
                c = abs(float(optimize.fsolve(self.objective, c, [epsc[i], epsc_prime, asinit], factor=100, xtol=1e-4)))
                # Stop analysis if RunTimeWarning is caught (i.e. no convergence)
                if math.isnan(self.mi):
                    # Check if target moment was reached (it not then analysis stopped prematurely due to bad guess)
                    if max(m) < self.m_target:
                        # Rerun with different c
                        c = 0.03
                        c = abs(float(optimize.fsolve(self.objective, c, [epsc[i], epsc_prime, asinit], factor=100,
                                                      xtol=1e-4)))
                    else:
                        # Check if c initial should be modified, as the analysis stopped prematurely
                        if m[-2] / m[-1] < 0.9:
                            c = 0.02
                            c = abs(float(optimize.fsolve(self.objective, c, [epsc[i], epsc_prime, asinit], factor=100,
                                                          xtol=1e-4)))
                        else:
                            if m[-2] / m[-1] < 0.9:
                                c = 0.1
                                c = abs(
                                    float(optimize.fsolve(self.objective, c, [epsc[i], epsc_prime, asinit], factor=100,
                                                          xtol=1e-4)))
                            else:
                                break

                # Sometimes even though the M-Phi is obtained, the solution is still converging, so we need to stop it
                # to avoid inverse slopes of Curvature
                if phi.size <= 1:
                    if phi[-1] > self.phii:
                        break

                # tensile reinforcement strains
                eps_tensile = np.append(eps_tensile, self.epss)
                # tensile reinforcement stresses
                sigmat = np.append(sigmat, self.fst)
                # bending moment capacity
                m = np.append(m, self.mi)
                # curvature
                phi = np.append(phi, self.phii)

                # Stop analysis if bottom reinforcement has ruptured
                # TODO, don't know why RESPONSE stops at half of strain,uk
                if self.epss >= self.EPSUK / 2:
                    break

            yield_index = getIndex(self.fy, sigmat)

            # Removing None arguments
            if self.k_hard == 1.:
                m = m[~np.isnan(m)]
                phi = phi[~np.isnan(phi)]
            else:
                idx = min(np.argwhere(np.isnan(m))[0][0], np.argwhere(np.isnan(phi))[0][0])
                m = m[:idx]
                phi = phi[:idx]
            idx_max = -1
            m_max = m[idx_max]
            my_first = m[yield_index]
            phiy_first = phi[yield_index]
            m = np.array(m)
            rpeak = m_max / my_first
            ei_cracked = my_first / phiy_first
            ei_cracked = ei_cracked / (young_modulus_rc * self.b * self.h ** 3 / 12 * 1000)

        # Nominal yield curvature
        phi_yield_nom = self.m_target * phiy_first / my_first

        # Curvature ductility (hardening ductility)
        mu_phi = phi[-1] / phi_yield_nom

        # todo, add accounting for residual strength near here
        # Softening slope (todo, Shear design to be added near here)
        nu = abs(self.p) / area / self.fc_prime / 1000
        ro_sh = self.TRANSVERSE_LEGS * np.pi * self.TRANSVERSE_DIAMETER ** 2 / 4 / self.TRANSVERSE_SPACING / self.b
        A_sh = self.TRANSVERSE_LEGS * np.pi * self.TRANSVERSE_DIAMETER ** 2 / 4

        if math.isnan(self.epss):
            self.epss = eps_tensile[-1]

        phi_critical, m_critical, lp = self.get_softening_slope(rebar_area=asinit, curvature_yield=phiy_first,
                                                                curvature_ductility=mu_phi, axial_load_ratio=nu,
                                                                transverse_steel_ratio=ro_sh)

        # M and Curvature for the lumped hinge model fitting / calibration if need be
        m_critical = 1e-9 if m_critical == 0 else m_critical
        m_model = np.array([1e-9, self.m_target, max(m), m_critical])

        # Force phi_critical equal to max phi, if not exceeding (for columns), it is unlikely to have this situation
        phi_critical = max(phi[-1], phi_critical)
        # TODO, check into phi_critical and get rid of the following if statement
        if phi_critical == max(phi):
            phi_critical = max(phi) * 1.01

        phi_model = np.array([1e-9, phi_yield_nom, max(phi), phi_critical])

        # Identifying fracturing point
        m = np.append(m, m_critical)
        phi = np.append(phi, phi_critical)
        fracturing_ductility = phi_critical / phi_yield_nom

        # Exporting the results
        # The values are relative to the point of yield definition (herein to first yield)
        # If columns are designed, full reinforcement is recorded, if beams, then the direction of interest is recorded
        As_factor = 1. if self.AsTotal is None else self.distAs[0]

        data = {'curvature': phi, 'moment': m, 'curvature_ductility': mu_phi, 'peak/yield ratio': rpeak,
                'reinforcement': asinit * As_factor, 'cracked EI': ei_cracked, 'first_yield_moment': my_first,
                'first_yield_curvature': phiy_first, 'phi_critical': phi_critical,
                'fracturing_ductility': fracturing_ductility, "lp": lp, "cover": self.d, "A_sh": A_sh,
                "spacing": self.TRANSVERSE_SPACING, "b": self.b, "h": self.h}
        reinforcement = {"Strain": eps_tensile, "Stress": sigmat}
        concrete = {"Strain": epsc, "Stress": sigma_c}
        MPhi_idealization = {"phi": phi_model, "m": m_model}

        # Hysteretic behaviour of all structural elements for model creation in OpenSees (M-curvature)
        # Assuming 50% of gross cross-section (the actual M-phi is calculated without the necessity of defining fstiff)
        # Those are useless, TO BE REMOVED
        curv_yield = self.m_target / young_modulus_rc / 1000 / inertia / self.fstiff
        curv_ult = mu_phi * phiy_first
        model = {"yield": {"curvature": curv_yield, "moment": self.m_target},
                 "ultimate": {"curvature": curv_ult, "moment": m_max},
                 "fracturing": {"curvature": phi_critical, "moment": 0}}

        return data, reinforcement, concrete, model, MPhi_idealization


if __name__ == '__main__':
    """
    --- Info on the input data:
    b                       section width [m]
    h                       section height [m]
    m_target                target Moment demand [kNm]
    nlayers                 number of reinforcement layers
    p                       analysis axial force, negative=compressive [kN]
    d                       reinforcement cover [m]
    fc_prime                concrete strength [MPa]
    fy                      reinforcement yield strength [MPa]
    young_modulus_s         reinforcement elastic modulus [MPa]
    """
    # Section properties
    b = 0.45
    h = 0.45
    Mtarget = 147.53
    N = 695.9
    cover = 0.03
    nlayers = 1

    mphi = MomentCurvatureRC(b, h, Mtarget, p=N, nlayers=nlayers, d=cover, young_mod_s=200000.,
                             k_hard=1, soft_method="Collins")

    m = mphi.get_mphi()
    print(m[3])

