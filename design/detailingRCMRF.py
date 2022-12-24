"""
defines detailing conditions (code-based) for element design
The detailing phase comes as a phase before an iteration where the SPO curve needs to be updated
"""
import numpy as np
from scipy import optimize
import pandas as pd

from design.momentcurvaturerc import MomentCurvatureRC
from design.plasticity import Plasticity


class DetailingRCMRF:
    def __init__(self, demands, nst, nbays, fy, fc, bay_widths, heights, n_seismic, mi, dy, sections,
                 rebar_cover=0.04, ductility_class="DCM", young_mod_s=200e3, k_hard=1.0, est_ductilities=True,
                 direction=0):
        """
        initializes detailing phase
        :param demands: dict                Demands on structural elements
        :param nst: int                     Number of stories
        :param nbays: int                   Number of bays
        :param fy: float                    Reinforcement yield strength
        :param fc: float                    Concrete compressive strength
        :param bay_widths: list             Bay widths
        :param heights: list                Storey heights
        :param n_seismic: int               Number of seismic frames
        :param mi: list                     Lumped storey masses
        :param dy: float                    System yield displacement in m
        :param sections: DataFrame          Cross-sections of elements of the solution
        :param rebar_cover: float           Reinforcement cover in m
        :param ductility_class: str         Ductility class (DCM or DCH, following Eurocode 8 recommendations)
        :param young_mod_s: float           Young modulus of reinforcement
        :param k_hard: float                Hardening slope of reinforcement (i.e. fu/fy)
        :param est_ductilities: bool        Whether to estimate global ductilities
        :param direction: int               Directivity of elements (X or Y), for 3D approach
        """
        self.demands = demands
        self.nst = nst
        self.nbays = nbays
        self.fy = fy
        self.young_mod_s = young_mod_s
        self.fc = fc
        self.bay_widths = bay_widths
        self.heights = heights
        self.n_seismic = n_seismic
        self.mi = mi
        self.dy = dy
        self.sections = sections
        self.rebar_cover = rebar_cover
        self.ductility_class = ductility_class
        # Reinforcement characteristic yield strength in MPa
        self.FYK = 500.
        self.k_hard = k_hard
        # if no warning, then no iterations are necessary for 4a
        self.WARNING_MAX = False
        self.WARNING_MIN = False
        # Warning for each element, if local max reinforcement ratio limit is not met
        self.WARN_ELE_MAX = False
        # Warning for each element, if local min reinforcement ratio limit is not met
        self.WARN_ELE_MIN = False
        # Estimate ductilities
        self.est_ductilities = est_ductilities
        # Direction of seismic action
        self.direction = direction

    def capacity_design(self, Mbi, Mci):
        """
        Applies capacity design strong column - weak beam concept.
        Assumption: Even though the peak moment might not occur with peak compressive or tension loads, the procedure
        still applies M+N and M-N with the peak moment, as the peak moment capacity is sometimes dictated by capacity
        design procedures.
        :param Mbi: ndarray                 Moment demands on beams
        :param Mci: ndarray                 Moment demands on columns
        :return: ndarray                    Beam and column moment demands
        """
        Myc = Mci.copy()
        Myb = Mbi.copy()

        for bay in range(self.nbays, -1, -1):
            # No capacity design for the top storey
            for st in range(self.nst - 2, -1, -1):
                if bay == 0:
                    diff = Myb[st][bay] * 1.3 - (Myc[st][bay] + Myc[st + 1][bay])
                    if diff > 0.0:
                        Myc[st][bay] += diff / 2
                        Myc[st + 1][bay] += diff / 2
                elif bay == self.nbays:
                    diff = Myb[st][bay - 1] * 1.3 - (Myc[st][bay] + Myc[st + 1][bay])
                    if diff > 0.0:
                        Myc[st][bay] += diff / 2
                        Myc[st + 1][bay] += diff / 2
                else:
                    diff = (Myb[st][bay - 1] + Myb[st][bay]) * 1.3 - (Myc[st][bay] + Myc[st + 1][bay])
                    if diff > 0.0:
                        Myc[st][bay] += diff / 2
                        Myc[st + 1][bay] += diff / 2
        return Myc

    def ensure_symmetry(self, option="max"):
        """
        Ensures symmetric strength along the width of the frame
        :param option: str                  Technique of reading demands ('max', 'mean', 'min')
        :return: ndarray                    Internal force demands of structural elements, i.e. M of beams,
                                            M and N of columns
        """
        # Beams
        MbiPos = self.demands["Beams"]["M"]["Pos"]
        MbiNeg = self.demands["Beams"]["M"]["Neg"]
        Mci = self.demands["Columns"]["M"]
        Nci = self.demands["Columns"]["N"]
        NciNeg = self.demands["Columns"]["N"].copy()
        if self.nbays <= 2:
            for st in range(self.nst):
                if option == "max":
                    MbiPos[st][0] = MbiPos[st][self.nbays - 1] = np.max((MbiPos[st][0], MbiPos[st][self.nbays - 1]))
                    MbiNeg[st][0] = MbiNeg[st][self.nbays - 1] = np.max((MbiNeg[st][0], MbiNeg[st][self.nbays - 1]))
                    Mci[st][0] = Mci[st][self.nbays] = np.max((Mci[st][0], Mci[st][self.nbays]))
                    Nci[st][0] = Nci[st][self.nbays] = np.max((Nci[st][0], Nci[st][self.nbays]))
                    NciNeg[st][0] = NciNeg[st][self.nbays] = min(np.min((NciNeg[st][0], NciNeg[st][self.nbays])), 0.0)
                elif option == "mean":
                    MbiPos[st][0] = MbiPos[st][self.nbays - 1] = np.mean((MbiPos[st][0], MbiPos[st][self.nbays - 1]))
                    MbiNeg[st][0] = MbiNeg[st][self.nbays - 1] = np.mean((MbiNeg[st][0], MbiNeg[st][self.nbays - 1]))
                    Mci[st][0] = Mci[st][self.nbays] = np.mean((Mci[st][0], Mci[st][self.nbays]))
                    Nci[st][0] = Nci[st][self.nbays] = np.mean((Nci[st][0], Nci[st][self.nbays]))
                    NciNeg[st][0] = NciNeg[st][self.nbays] = min(np.min((NciNeg[st][0], NciNeg[st][self.nbays])), 0.0)
                elif option == "min":
                    MbiPos[st][0] = MbiPos[st][self.nbays - 1] = np.min((MbiPos[st][0], MbiPos[st][self.nbays - 1]))
                    MbiNeg[st][0] = MbiNeg[st][self.nbays - 1] = np.min((MbiNeg[st][0], MbiNeg[st][self.nbays - 1]))
                    Mci[st][0] = Mci[st][self.nbays] = np.min((Mci[st][0], Mci[st][self.nbays]))
                    Nci[st][0] = Nci[st][self.nbays] = np.min((Nci[st][0], Nci[st][self.nbays]))
                    NciNeg[st][0] = NciNeg[st][self.nbays] = min(np.min((NciNeg[st][0], NciNeg[st][self.nbays])), 0.0)
                else:
                    raise ValueError("[EXCEPTION] Wrong option for ensuring symmetry, must be max, mean or min")
        else:
            for st in range(self.nst):

                if option == "max":
                    MbiPos[st][0] = MbiPos[st][self.nbays - 1] = np.max((MbiPos[st][0], MbiPos[st][self.nbays - 1]))
                    MbiPos[st][1:self.nbays - 1] = np.max(MbiPos[st][1:self.nbays - 1])
                    MbiNeg[st][0] = MbiNeg[st][self.nbays - 1] = np.max((MbiNeg[st][0], MbiNeg[st][self.nbays - 1]))
                    MbiNeg[st][1:self.nbays - 1] = np.max(MbiNeg[st][1:self.nbays - 1])
                    Mci[st][0] = Mci[st][self.nbays] = np.max((Mci[st][0], Mci[st][self.nbays]))
                    for bay in range(1, int((self.nbays - 1) / 2) + 1):
                        Mci[st][bay] = Mci[st][-bay - 1] = np.max((Mci[st][bay], Mci[st][-bay - 1]))
                    Nci[st][0] = Nci[st][self.nbays] = np.max((Nci[st][0], Nci[st][self.nbays]))
                    NciNeg[st][0] = NciNeg[st][self.nbays] = min(np.min((NciNeg[st][0], NciNeg[st][self.nbays])), 0.0)
                    for bay in range(1, int((self.nbays - 1) / 2) + 1):
                        Nci[st][bay] = Nci[st][-bay - 1] = np.max((Nci[st][bay], Nci[st][-bay - 1]))
                        NciNeg[st][bay] = NciNeg[st][-bay - 1] = min(np.min((NciNeg[st][bay], NciNeg[st][-bay - 1])),
                                                                     0.0)

                elif option == "mean":
                    MbiPos[st][0] = MbiPos[st][self.nbays - 1] = np.mean((MbiPos[st][0], MbiPos[st][self.nbays - 1]))
                    MbiPos[st][1:self.nbays - 1] = np.mean(MbiPos[st][1:self.nbays - 1])
                    MbiNeg[st][0] = MbiNeg[st][self.nbays - 1] = np.mean((MbiNeg[st][0], MbiNeg[st][self.nbays - 1]))
                    MbiNeg[st][1:self.nbays - 1] = np.mean(MbiNeg[st][1:self.nbays - 1])
                    Mci[st][0] = Mci[st][self.nbays] = np.mean((Mci[st][0], Mci[st][self.nbays]))
                    for bay in range(1, int((self.nbays - 1) / 2) + 1):
                        Mci[st][bay] = Mci[st][-bay - 1] = np.mean((Mci[st][bay], Mci[st][-bay - 1]))
                    Nci[st][0] = Nci[st][self.nbays] = np.mean((Nci[st][0], Nci[st][self.nbays]))
                    NciNeg[st][0] = NciNeg[st][self.nbays] = min(np.min((NciNeg[st][0], NciNeg[st][self.nbays])), 0.0)
                    for bay in range(1, int((self.nbays - 1) / 2) + 1):
                        Nci[st][bay] = Nci[st][-bay - 1] = np.mean((Nci[st][bay], Nci[st][-bay - 1]))
                        NciNeg[st][bay] = NciNeg[st][-bay - 1] = min(np.min((NciNeg[st][bay], NciNeg[st][-bay - 1])),
                                                                     0.0)

                elif option == "min":
                    MbiPos[st][0] = MbiPos[st][self.nbays - 1] = np.min((MbiPos[st][0], MbiPos[st][self.nbays - 1]))
                    MbiPos[st][1:self.nbays - 1] = np.min(MbiPos[st][1:self.nbays - 1])
                    MbiNeg[st][0] = MbiNeg[st][self.nbays - 1] = np.min((MbiNeg[st][0], MbiNeg[st][self.nbays - 1]))
                    MbiNeg[st][1:self.nbays - 1] = np.min(MbiNeg[st][1:self.nbays - 1])
                    Mci[st][0] = Mci[st][self.nbays] = np.min((Mci[st][0], Mci[st][self.nbays]))
                    for bay in range(1, int((self.nbays - 1) / 2) + 1):
                        Mci[st][bay] = Mci[st][-bay - 1] = np.min((Mci[st][bay], Mci[st][-bay - 1]))
                    Nci[st][0] = Nci[st][self.nbays] = np.min((Nci[st][0], Nci[st][self.nbays]))
                    NciNeg[st][0] = NciNeg[st][self.nbays] = min(np.min((NciNeg[st][0], NciNeg[st][self.nbays])), 0.0)
                    for bay in range(1, int((self.nbays - 1) / 2) + 1):
                        Nci[st][bay] = Nci[st][-bay - 1] = np.min((Nci[st][bay], Nci[st][-bay - 1]))
                        NciNeg[st][bay] = NciNeg[st][-bay - 1] = min(np.min((NciNeg[st][bay], NciNeg[st][-bay - 1])),
                                                                     0.0)
                else:
                    raise ValueError("[EXCEPTION] Wrong option for ensuring symmetry, must be max, mean or min")
        return MbiPos, MbiNeg, Mci, Nci, NciNeg

    def ensure_local_ductility(self, b, h, reinforcement, relation, st, bay, eletype, oppReinf=None, pflag=True):
        """
        Local ductility checks according to Eurocode 8
        :param b: float                             Width of element
        :param h: float                             Height of element
        :param reinforcement: float                 Total reinforcement area
        :param relation: class                      Object created based on M-phi relationship
        :param st: int                              Storey level
        :param bay: int                             Bay level
        :param eletype: str                         Element type, beam or column
        :param oppReinf: float                      Reinforcement of opposite direction
        :param pflag: bool                          Print out warnings
        :return: dict                               M-phi outputs to be stored
        """
        # Behaviour factor, for frame systems assuming regularity in elevation
        # Assuming multi-storey, multi-bay frames
        if self.ductility_class == "DCM":
            q = 3.9
        elif self.ductility_class == "DCH":
            q = 5.85
        else:
            raise ValueError("[EXCEPTION] Wrong type of ductility class, must be DCM or DCH!")

        # Design compressive strength
        fcd = 0.8 / 1.5 * self.fc
        # Tensile strength
        fctm = 0.3 * self.fc ** (2 / 3)

        # Reinforcement ratio
        ro_prime = reinforcement / (b * (h - self.rebar_cover))
        req_duct = q * 2 - 1
        yield_strain = self.fy / self.young_mod_s

        # Reinforcement ratio limits to meet local ductility conditions
        if eletype == "Beam":
            ro_max = ro_prime + 0.0018 * fcd / self.fy / req_duct / yield_strain
            ro_min = 0.5 * fctm / self.FYK
        elif eletype == "Column":
            ro_min = 0.01
            ro_max = 0.04
        else:
            raise ValueError("[EXCEPTION] Wrong type of element for local ductility verifications!")

        # Verifications
        if ro_min > ro_prime:
            ro_prime = ro_min
            if eletype == "Beam":
                # self.WARN_ELE_MIN = True
                # cover = relation.d
                # while self.WARN_ELE_MIN and cover < 0.04 - 0.0005:
                    # Increase reinforcement cover, which will trigger requirement of more reinforcement
                #     cover += 0.005
                #     data = relation.get_mphi(cover=cover)
                #     reinforcement = data[0]["reinforcement"]
                #     ro_prime = reinforcement / (b * (h - cover))
                #
                #     if ro_min > ro_prime:
                #         self.WARN_ELE_MIN = True
                #     else:
                #         self.WARN_ELE_MIN = False
                # data = relation.get_mphi(cover=cover)

                rebar = (b * (h - self.rebar_cover)) * ro_prime + oppReinf
                m_target = relation.get_mphi(check_reinforcement=True, reinf_test=rebar,
                                             reinforcements=[(b * (h - self.rebar_cover)) * ro_prime, oppReinf])
                data = relation.get_mphi(m_target=m_target,
                                         reinforcements=[(b * (h - self.rebar_cover)) * ro_prime, oppReinf])
            else:
                # while self.WARN_ELE_MIN:
                #     # Increase reinforcement cover, which will trigger requirement of more reinforcement
                #     cover = relation.d + 0.05
                #     data = relation.get_mphi(cover=cover)

                rebar = (b * (h - self.rebar_cover)) * ro_prime
                m_target = relation.get_mphi(check_reinforcement=True, reinf_test=rebar)
                data = relation.get_mphi(m_target=m_target)

            self.WARN_ELE_MAX = False
            self.WARN_ELE_MIN = True
            self.WARNING_MIN = True
            return data

        elif ro_max < ro_prime:
            if pflag:
                print(f"[WARNING] Cross-section of {eletype} element at storey {st} and bay {bay} should be increased! "
                      f"ratio: {ro_prime * 100:.2f}%")
                self.WARN_ELE_MAX = True
                self.WARN_ELE_MIN = False
                self.WARNING_MAX = True
            return None

        else:
            self.WARN_ELE_MIN = False
            self.WARN_ELE_MAX = False
            return None

    def compute_m_based_on_reinforcement(self, b, h, d, As):
        eta = min(1., 1 - (self.fc - 50) / 200)
        alpha_cc = 0.85
        # Partial factor for concrete at ULS assuming persistent and transient design situations
        gamma_uls = 1.5
        fcd = alpha_cc * self.fc / gamma_uls * 1000
        fy = self.fy * 1000
        return fy * As * (h - d) * (1 - fy * As / (2 * (h - d) * fcd * eta * b))

    def get_rebar_distribution(self, b, h, d, mpos, mneg):
        """
        Gets initial rebar distribution based on Eurocode design procedures.
        Gives an initial guess on the proportions and values of reinforcement in both sides of the beams.
        :param b: float                         Widths of element
        :param h: float                         Height of element
        :param d: float                         Effective height of element
        :param mpos: float                      Positive moment demand
        :param mneg: float                      Negative moment demand
        :return: float, list                    Total reinforcement area in m2 and relative distributions of reinforcement
        """
        eta = min(1., 1 - (self.fc - 50) / 200)
        alpha_cc = 0.85
        # Partial factor for concrete at ULS assuming persistent and transient design situations
        gamma_uls = 1.5
        fcd = alpha_cc * self.fc / gamma_uls * 1000
        fy = self.fy * 1000

        def get_As(As, moment):
            return fy * As * (h - d) * (1 - fy * As / (2 * (h - d) * fcd * eta * b)) - moment

        # Initial guess for the solver
        As = 0.002
        AsPos = float(optimize.fsolve(get_As, As, mpos, factor=0.1))
        AsNeg = float(optimize.fsolve(get_As, As, mneg, factor=0.1))
        AsTotal = AsPos + AsNeg
        distributions = [AsPos / AsTotal, AsNeg / AsTotal]
        return AsTotal, distributions

    def design_gravity(self):
        """
        designs gravity elements using demands from ELFM and optimal solution, uses moment_curvature_rc
        :return: dict                           Designed element details from the moment-curvature relationship
        """
        # Demands on all structural elements
        demands = self.demands
        # Get only the peak values for moments and moments for beams
        # Assumption: all beams independent of direction within a storey are of the same section
        b_pos = np.zeros((self.nst, 2))
        b_neg = np.zeros((self.nst, 2))
        for st in range(self.nst):
            cnt = 0
            for eletype in demands["x"]:
                if eletype.startswith("Beam"):
                    # Positive and negative sides
                    for s in demands["x"][eletype]["M"]:
                        tempx = np.max(np.abs(demands["x"][eletype]["M"][s][st, :, :]))
                        tempy = np.max(np.abs(demands["y"][eletype]["M"][s][st, :, :]))
                        if s == "Pos":
                            b_pos[st][cnt] = max(tempx, tempy)
                        else:
                            b_neg[st][cnt] = max(tempx, tempy)

                    cnt += 1

        # Demands on central/gravity columns. Get only the peak demands
        # Assumption: Columns within a storey are the same within the storey
        c = np.zeros((self.nst, 2))
        for st in range(self.nst):
            # Columns
            col_x = np.max(np.abs(demands["x"]["Columns"]["M"][st, :, :]))
            col_y = np.max(np.abs(demands["y"]["Columns"]["M"][st, :, :]))
            if col_x > col_y:
                c[st][0] = col_x
                idxs = list(zip(*np.where(np.abs(demands["x"]["Columns"]["M"][st, :, :]) == c[st][0])))[0]
                c[st][1] = np.abs(demands["x"]["Columns"]["N"][st, idxs[0], idxs[1]])
            else:
                c[st][0] = col_y
                idxs = list(zip(*np.where(np.abs(demands["y"]["Columns"]["M"][st, :, :]) == c[st][0])))[0]
                c[st][1] = np.abs(demands["y"]["Columns"]["N"][st, idxs[0], idxs[1]])

        # Maximum moment demands for beams at each storey level
        # Different beams along each of the directions
        beam_demands_pos = b_pos
        beam_demands_neg = b_neg

        # Column demand at each storey level. Moment; Axial load
        column_demands = c

        # Initialize hinge models
        model = {"Beams": {"Pos": {}, "Neg": {}}, "Columns": {}}
        hinge_models = {"Beams": {"Pos": {}, "Neg": {}}, "Columns": {}}
        warnings = {"MAX": {"Beams": {"Pos": {}, "Neg": {}}, "Columns": {}},
                    "MIN": {"Beams": {"Pos": {}, "Neg": {}}, "Columns": {}}}

        # Design of beams
        for st in range(self.nst):
            # Initialize for each storey
            model["Beams"]["Pos"][f"S{st + 1}"] = {"x": {}, "y": {}}
            model["Beams"]["Neg"][f"S{st + 1}"] = {"x": {}, "y": {}}
            hinge_models["Beams"]["Pos"][f"S{st + 1}"] = {"x": {}, "y": {}}
            hinge_models["Beams"]["Neg"][f"S{st + 1}"] = {"x": {}, "y": {}}

            # Along x direction for space systems
            m_target_pos = beam_demands_pos[st][0]
            m_target_neg = beam_demands_neg[st][0]
            b = self.sections[f"bx{st + 1}"]
            h = self.sections[f"hx{st + 1}"]
            model_pos, hinge_pos, model_neg, hinge_neg, w = self.get_details(b, h, m_target_pos, m_target_neg)
            model["Beams"]["Pos"][f"S{st + 1}"]["x"] = model_pos
            model["Beams"]["Neg"][f"S{st + 1}"]["x"] = model_neg
            hinge_models["Beams"]["Pos"][f"S{st + 1}"]["x"] = hinge_pos
            hinge_models["Beams"]["Neg"][f"S{st + 1}"]["x"] = hinge_neg
            # Add the warnings
            warnings["MAX"]["Beams"]["Pos"][f"S{st + 1}x"] = w["pos"][0]
            warnings["MIN"]["Beams"]["Pos"][f"S{st + 1}x"] = w["pos"][1]
            warnings["MAX"]["Beams"]["Neg"][f"S{st + 1}x"] = w["neg"][0]
            warnings["MIN"]["Beams"]["Neg"][f"S{st + 1}x"] = w["neg"][1]

            # Detail elements along Y direction for space systems
            m_target_pos = beam_demands_pos[st][1]
            m_target_neg = beam_demands_neg[st][1]
            b = self.sections[f"by{st + 1}"]
            h = self.sections[f"hy{st + 1}"]
            model_pos, hinge_pos, model_neg, hinge_neg, w = self.get_details(b, h, m_target_pos, m_target_neg)
            model["Beams"]["Pos"][f"S{st + 1}"]["y"] = model_pos
            model["Beams"]["Neg"][f"S{st + 1}"]["y"] = model_neg
            hinge_models["Beams"]["Pos"][f"S{st + 1}"]["y"] = hinge_pos
            hinge_models["Beams"]["Neg"][f"S{st + 1}"]["y"] = hinge_neg
            # Add the warnings
            warnings["MAX"]["Beams"]["Pos"][f"S{st + 1}y"] = w["pos"][0]
            warnings["MIN"]["Beams"]["Pos"][f"S{st + 1}y"] = w["pos"][1]
            warnings["MAX"]["Beams"]["Neg"][f"S{st + 1}y"] = w["neg"][0]
            warnings["MIN"]["Beams"]["Neg"][f"S{st + 1}y"] = w["neg"][1]

        # Design of columns
        for st in range(self.nst):
            b = h = self.sections[f"hi{st + 1}"]
            # Design bending moment
            m_target_pos = column_demands[st][0]
            # Design compressive internal axial force
            nc_design = column_demands[st][1]

            # Number of reinforcement layers based on section height (may be adjusted manually)
            nlayers = 0 if h <= 0.3 else 1 if (0.3 < h <= 0.55) else 2
            # Assuming contraflexure at 0.6 of height
            z = 0.6 * self.heights[st]

            mphi = MomentCurvatureRC(b, h, m_target_pos, length=z, p=-nc_design, nlayers=nlayers, d=self.rebar_cover,
                                     young_mod_s=self.young_mod_s, k_hard=self.k_hard, soft_method="Collins",
                                     fy=self.fy, fc_prime=self.fc)

            model["Columns"][f"S{st + 1}"] = mphi.get_mphi()
            hinge_models["Columns"][f"S{st + 1}"] = model["Columns"][f"S{st + 1}"][4]

            '''Local ductility requirement checks (following Eurocode 8 recommendations)'''
            model_temp = self.ensure_local_ductility(b, h, model["Columns"][f"S{st + 1}"][0]["reinforcement"], mphi,
                                                     None, None, eletype="Column", pflag=False)
            # Any warnings
            warnings["MAX"]["Columns"][f"S{st + 1}"] = self.WARN_ELE_MAX
            warnings["MIN"]["Columns"][f"S{st + 1}"] = self.WARN_ELE_MIN

            if model_temp is not None:
                model["Columns"][f"S{st + 1}"] = model_temp
                hinge_models["Columns"][f"S{st + 1}"] = model_temp[4]

        # Get hinge model information in DataFrame
        columns = ["Element", "Storey", "Direction", "b", "h", "coverNeg", "coverPos", "lp", "phi1Neg", "phi2Neg",
                   "phi3Neg", "m1Neg", "m2Neg", "m3Neg", "phi1", "phi2", "phi3", "m1", "m2", "m3"]

        numericCols = ["Storey", "b", "h", "Direction", "coverNeg", "coverPos", "lp", "phi1Neg", "phi2Neg", "phi3Neg",
                       "m1Neg", "m2Neg", "m3Neg", "phi1", "phi2", "phi3", "m1", "m2", "m3"]

        hinge_models = pd.DataFrame(columns=columns)
        for ele in model:
            m = model[ele]
            if ele.lower() == "beams":
                for j in m["Pos"]:
                    lp = m["Pos"][j]["x"][0]["lp"]
                    st = int(j[1])
                    # Along X direction
                    temp = np.array([ele[:-1], st, 0, model[ele]["Neg"][j]["x"][0]["b"],
                                     model[ele]["Neg"][j]["x"][0]["h"], model[ele]["Neg"][j]["x"][0]["cover"],
                                     model[ele]["Pos"][j]["x"][0]["cover"], lp])
                    phiNeg = model[ele]["Neg"][j]["x"][4]["phi"][1:]
                    mNeg = model[ele]["Neg"][j]["x"][4]["m"][1:]
                    phiPos = model[ele]["Pos"][j]["x"][4]["phi"][1:]
                    mPos = model[ele]["Pos"][j]["x"][4]["m"][1:]

                    data = np.concatenate((temp, phiNeg, mNeg, phiPos, mPos)).reshape(1, len(columns))
                    # Concatenate into the DataFrame
                    hinge_models = hinge_models.append(pd.DataFrame(data=data, columns=columns), ignore_index=True)

                    # Along Y direction
                    lp = m["Pos"][j]["y"][0]["lp"]
                    temp = np.array([ele[:-1], st, 1, model[ele]["Neg"][j]["y"][0]["b"],
                                     model[ele]["Neg"][j]["y"][0]["h"], model[ele]["Neg"][j]["y"][0]["cover"],
                                     model[ele]["Pos"][j]["y"][0]["cover"], lp])
                    phiNeg = model[ele]["Neg"][j]["y"][4]["phi"][1:]
                    mNeg = model[ele]["Neg"][j]["y"][4]["m"][1:]
                    phiPos = model[ele]["Pos"][j]["y"][4]["phi"][1:]
                    mPos = model[ele]["Pos"][j]["y"][4]["m"][1:]

                    data = np.concatenate((temp, phiNeg, mNeg, phiPos, mPos)).reshape(1, len(columns))
                    # Concatenate into the DataFrame
                    hinge_models = hinge_models.append(pd.DataFrame(data=data, columns=columns), ignore_index=True)

            else:
                for j in m:
                    lp = m[j][0]["lp"]
                    st = int(j[1])
                    temp = np.array([ele[:-1], st, 0, model[ele][j][0]["b"], model[ele][j][0]["h"],
                                     model[ele][j][0]["cover"], model[ele][j][0]["cover"], lp])
                    phiNeg = phiPos = model[ele][j][4]["phi"][1:]
                    mNeg = mPos = model[ele][j][4]["m"][1:]

                    data = np.concatenate((temp, phiNeg, mNeg, phiPos, mPos)).reshape(1, len(columns))

                    # Concatenate into the DataFrame
                    hinge_models = hinge_models.append(pd.DataFrame(data=data, columns=columns), ignore_index=True)

        for i in numericCols:
            hinge_models[i] = pd.to_numeric(hinge_models[i])

        return model, hinge_models, warnings

    def design_elements(self, modes=None):
        """
        designs elements using demands from ELFM and optimal solution, uses moment_curvature_rc
        :param modes: dict                      Periods and modal shapes obtained from modal analysis
        :return: dict                           Designed element details from the moment-curvature relationship
        """
        # Ensure symmetry of strength distribution along the widths of the frame
        ''' Assumptions: Beams are designed for both directions: positive and negative, neglecting axial loads
        Columns are designed considering M+N demands; No shear design is carried out'''
        mbiPos, mbiNeg, mci, nci, nciNeg = self.ensure_symmetry(option="max")

        # Initialize dictionaries for storing details, warnings and hinge models (hysteretic models) for OpenSees model
        data = {"Beams": {"Pos": {}, "Neg": {}}, "Columns": {}}
        warnings = {"MAX": {"Beams": {"Pos": {}, "Neg": {}}, "Columns": {}},
                    "MIN": {"Beams": {"Pos": {}, "Neg": {}}, "Columns": {}}}
        hinge_models = {"Beams": {"Pos": {}, "Neg": {}}, "Columns": {}}

        # Design of beams
        for st in range(self.nst):
            if self.nbays > 2:
                for bay in range(int(round(self.nbays / 2, 0))):
                    # Design bending moment
                    # Note: Negative = bottom, positive = top
                    m_target_pos = mbiPos[st][bay]
                    m_target_neg = mbiNeg[st][bay]

                    # Cross-section dimensions
                    b = self.sections[f"b{st + 1}"]
                    h = self.sections[f"h{st + 1}"]

                    # Initial guess on the distribution and values of the reinforcements
                    AsTotal, distributions = self.get_rebar_distribution(b, h, self.rebar_cover, m_target_pos,
                                                                         m_target_neg)

                    # TODO, modify so that Negative direction is run with the knowledge of AsPos and seeks only AsNeg
                    # Perform moment-curvature analysis, Positive direction
                    mphiPos = MomentCurvatureRC(b, h, m_target_pos, d=self.rebar_cover, young_mod_s=self.young_mod_s,
                                                k_hard=self.k_hard, AsTotal=AsTotal, distAs=distributions,
                                                fy=self.fy, fc_prime=self.fc)
                    data["Beams"]["Pos"][f"S{st + 1}B{bay + 1}"] = mphiPos.get_mphi()

                    # Negative direction
                    mphiNeg = MomentCurvatureRC(b, h, m_target_neg, d=self.rebar_cover, young_mod_s=self.young_mod_s,
                                                k_hard=self.k_hard, AsTotal=AsTotal, distAs=distributions[::-1],
                                                fy=self.fy, fc_prime=self.fc)
                    data["Beams"]["Neg"][f"S{st + 1}B{bay + 1}"] = mphiNeg.get_mphi()

                    # Hinge models
                    hinge_models["Beams"]["Pos"][f"S{st+1}B{bay+1}"] = data["Beams"]["Pos"][f"S{st+1}B{bay+1}"][4]
                    hinge_models["Beams"]["Neg"][f"S{st+1}B{bay+1}"] = data["Beams"]["Neg"][f"S{st+1}B{bay+1}"][4]

                    '''Local ductility requirement checks (following Eurocode 8 recommendations)'''
                    # Positive direction
                    d_temp = self.ensure_local_ductility(b, h, data["Beams"]["Pos"][f"S{st + 1}B{bay + 1}"][0][
                        "reinforcement"],  mphiPos, st + 1, bay + 1, eletype="Beam",
                                                         oppReinf=data["Beams"]["Neg"][f"S{st + 1}B{bay + 1}"][0][
                                                             "reinforcement"])

                    warnings["MAX"]["Beams"]["Pos"][f"S{st + 1}B{bay + 1}"] = self.WARN_ELE_MAX
                    warnings["MIN"]["Beams"]["Pos"][f"S{st + 1}B{bay + 1}"] = self.WARN_ELE_MIN
                    if d_temp is not None:
                        data["Beams"]["Pos"][f"S{st + 1}B{bay + 1}"] = d_temp
                        hinge_models["Beams"]["Pos"][f"S{st+1}B{bay+1}"] = d_temp[4]

                    # Negative direction
                    d_temp = self.ensure_local_ductility(b, h, data["Beams"]["Neg"][f"S{st + 1}B{bay + 1}"][0][
                        "reinforcement"], mphiNeg, st + 1, bay + 1, eletype="Beam",
                                                         oppReinf=data["Beams"]["Pos"][f"S{st + 1}B{bay + 1}"][0][
                                                             "reinforcement"])

                    # TODO, once local ductility is ensured, M-phi relationship might change, also after mphiNeg, pos reinforcement might change
                    # So, ideally it should go back and forth to correct the reinforcements, however, no iterations are done there
                    warnings["MAX"]["Beams"]["Neg"][f"S{st + 1}B{bay + 1}"] = self.WARN_ELE_MAX
                    warnings["MIN"]["Beams"]["Neg"][f"S{st + 1}B{bay + 1}"] = self.WARN_ELE_MIN
                    if d_temp is not None:
                        data["Beams"]["Neg"][f"S{st + 1}B{bay + 1}"] = d_temp
                        hinge_models["Beams"]["Neg"][f"S{st+1}B{bay+1}"] = d_temp[4]

                    # Update the beam demands for capacity design of columns
                    mbiPos[st][bay] = data['Beams']['Pos'][f"S{st + 1}B{bay + 1}"][4]["m"][1]
                    mbiNeg[st][bay] = data['Beams']['Neg'][f"S{st + 1}B{bay + 1}"][4]["m"][1]

            else:
                # Design bending moment
                m_target_pos = mbiPos[st][0]
                m_target_neg = mbiNeg[st][0]

                # Cross-section dimensions
                b = self.sections[f"b{st + 1}"]
                h = self.sections[f"h{st + 1}"]

                # Initial guess on the distribution and values of the reinforcements
                AsTotal, distributions = self.get_rebar_distribution(b, h, self.rebar_cover, m_target_pos, m_target_neg)

                # Perform moment-curvature analysis, Positive direction
                mphiPos = MomentCurvatureRC(b, h, m_target_pos, d=self.rebar_cover, young_mod_s=self.young_mod_s,
                                            k_hard=self.k_hard, AsTotal=AsTotal, distAs=distributions,
                                            fy=self.fy, fc_prime=self.fc)
                data["Beams"]["Pos"][f"S{st + 1}B{1}"] = mphiPos.get_mphi()
                # Negative direction
                mphiNeg = MomentCurvatureRC(b, h, m_target_neg, d=self.rebar_cover, young_mod_s=self.young_mod_s,
                                            k_hard=self.k_hard, AsTotal=AsTotal, distAs=distributions[::-1],
                                            fy=self.fy, fc_prime=self.fc)
                data["Beams"]["Neg"][f"S{st + 1}B{1}"] = mphiNeg.get_mphi()

                # Hinge models
                hinge_models["Beams"]["Pos"][f"S{st+1}B{1}"] = data["Beams"]["Pos"][f"S{st+1}B{1}"][4]
                hinge_models["Beams"]["Neg"][f"S{st+1}B{1}"] = data["Beams"]["Neg"][f"S{st+1}B{1}"][4]

                '''Local ductility requirement checks (following Eurocode 8 recommendations)'''
                # Positive direction
                d_temp = self.ensure_local_ductility(b, h, data["Beams"]["Pos"][f"S{st + 1}B{1}"][0]["reinforcement"],
                                                     mphiPos, st + 1, 1, eletype="Beam",
                                                     oppReinf=data["Beams"]["Neg"][f"S{st + 1}B{1}"][0][
                                                         "reinforcement"])
                warnings["MAX"]["Beams"]["Pos"][f"S{st + 1}B{1}"] = self.WARN_ELE_MAX
                warnings["MIN"]["Beams"]["Pos"][f"S{st + 1}B{1}"] = self.WARN_ELE_MIN
                if d_temp is not None:
                    data["Beams"]["Pos"][f"S{st + 1}B{1}"] = d_temp
                    hinge_models["Beams"]["Pos"][f"S{st+1}B{1}"] = d_temp[4]

                # Negative direction
                d_temp = self.ensure_local_ductility(b, h, data["Beams"]["Neg"][f"S{st + 1}B{1}"][0]["reinforcement"],
                                                     mphiNeg, st + 1, 1, eletype="Beam",
                                                     oppReinf=data["Beams"]["Pos"][f"S{st + 1}B{1}"][0][
                                                         "reinforcement"])
                warnings["MAX"]["Beams"]["Neg"][f"S{st + 1}B{1}"] = self.WARN_ELE_MAX
                warnings["MIN"]["Beams"]["Neg"][f"S{st + 1}B{1}"] = self.WARN_ELE_MIN
                if d_temp is not None:
                    data["Beams"]["Neg"][f"S{st + 1}B{1}"] = d_temp
                    hinge_models["Beams"]["Neg"][f"S{st+1}B{1}"] = d_temp[4]

                # Update the beam demands for capacity design of columns
                mbiPos[st][0] = data['Beams']['Pos'][f"S{st + 1}B{1}"][4]["m"][1]
                mbiNeg[st][0] = data['Beams']['Neg'][f"S{st + 1}B{1}"][4]["m"][1]

        ''' Follow the capacity design requirements, currently only Eurocode 8, but it may be adapted for other codes
        as well. Maximum absolute value of positive and negative moment demands on beams is used.'''
        # Max value to be used in capacity design
        mbi = np.maximum(mbiPos, mbiNeg)

        # Demands on columns
        myc = self.capacity_design(mbi, mci)

        # Design of columns
        for st in range(self.nst):
            for bay in range(int(np.ceil((self.nbays + 1) / 2))):
                if bay == 0:
                    b = h = self.sections[f"he{st + 1}"]
                else:
                    b = h = self.sections[f"hi{st + 1}"]
                # Design bending moment
                m_target = myc[st][bay]
                # Design compressive internal axial force
                nc_design = nci[st][bay]

                # Get tensile internal axial forces
                nc_design_neg = nciNeg[st][bay]

                # Number of reinforcement layers based on section height (may be adjusted manually)
                nlayers = 0 if h <= 0.3 else 1 if (0.3 < h <= 0.55) else 2
                # Assuming contraflexure at 0.6 of height
                # todo, may add better estimation of contraflexure point based on Muto's approach
                # todo, Collins softening method not working well with columns
                z = 0.6 * self.heights[st]

                mphi = MomentCurvatureRC(b, h, m_target, length=z, p=-nc_design, nlayers=nlayers, d=self.rebar_cover,
                                         young_mod_s=self.young_mod_s, k_hard=self.k_hard, soft_method="Collins",
                                         fy=self.fy, fc_prime=self.fc)

                temp = {"Pos": mphi.get_mphi()}

                if nc_design_neg < 0.0:
                    mphiNeg = MomentCurvatureRC(b, h, m_target, length=z, p=-nc_design_neg, nlayers=nlayers,
                                                d=self.rebar_cover, young_mod_s=self.young_mod_s, k_hard=self.k_hard,
                                                soft_method="Collins", fy=self.fy, fc_prime=self.fc)
                    temp["Neg"] = mphiNeg.get_mphi()
                    # Select the design requiring highest reinforcement
                    if temp["Neg"][0]["reinforcement"] > temp["Pos"][0]["reinforcement"]:
                        selection = temp["Neg"]
                        mphi = mphiNeg
                    else:
                        selection = temp["Pos"]
                else:
                    selection = temp["Pos"]

                data["Columns"][f"S{st + 1}B{bay + 1}"] = selection
                hinge_models["Columns"][f"S{st + 1}B{bay + 1}"] = selection[4]

                '''Local ductility requirement checks (following Eurocode 8 recommendations)'''
                d_temp = self.ensure_local_ductility(b, h, data["Columns"][f"S{st + 1}B{bay + 1}"][0]["reinforcement"],
                                                     mphi, st + 1, bay + 1, eletype="Column")
                warnings["MAX"]["Columns"][f"S{st + 1}B{bay + 1}"] = self.WARN_ELE_MAX
                warnings["MIN"]["Columns"][f"S{st + 1}B{bay + 1}"] = self.WARN_ELE_MIN
                if d_temp is not None:
                    data["Columns"][f"S{st + 1}B{bay + 1}"] = d_temp
                    hinge_models["Columns"][f"S{st + 1}B{bay + 1}"] = d_temp[4]

        # Old version, requires improvement
        if self.est_ductilities:
            mu_c, mu_f = self.estimate_ductilities(data, modes)
        else:
            mu_c = mu_f = None

        # Get hinge model information in DataFrame
        hinge_models = self.model_to_df(data)

        return data, hinge_models, mu_c, mu_f, warnings

    def get_details(self, b, h, m_pos, m_neg):
        """
        Gets details
        Note:   Update all detailing functions to use this one
                Currently run only for central elements (gravity)
        :param b: float                     Widths of element
        :param h: float                     Height of element
        :param m_pos: float                 Positive bending moment
        :param m_neg: float                 Negative bending moment
        :return: dictionaries               Details
        """
        AsTotal, distributions = self.get_rebar_distribution(b, h, self.rebar_cover, m_pos, m_neg)
        # Positive direction
        mphiPos = MomentCurvatureRC(b, h, m_pos, d=self.rebar_cover, young_mod_s=self.young_mod_s,
                                    k_hard=self.k_hard, AsTotal=AsTotal, distAs=distributions,
                                    fy=self.fy, fc_prime=self.fc)
        model_pos = mphiPos.get_mphi()

        # Negative direction
        mphiNeg = MomentCurvatureRC(b, h, m_neg, d=self.rebar_cover, young_mod_s=self.young_mod_s,
                                    k_hard=self.k_hard, AsTotal=AsTotal, distAs=distributions[::-1],
                                    fy=self.fy, fc_prime=self.fc)
        model_neg = mphiNeg.get_mphi()

        # Update the hinge models
        hinge_models_pos = model_pos
        hinge_models_neg = model_neg

        # Ensure local ductility in positive direction
        model_temp = self.ensure_local_ductility(b, h, model_pos[0]["reinforcement"], mphiPos, None, None,
                                                 eletype="Beam", oppReinf=model_neg[0]["reinforcement"], pflag=False)
        # Any warnings
        w = {"pos": [self.WARN_ELE_MAX, self.WARN_ELE_MIN]}

        if model_temp is not None:
            model_pos = model_temp
            hinge_models_pos = model_temp[4]

        # Ensure local ductility in negative direction
        model_temp = self.ensure_local_ductility(b, h, model_neg[0]["reinforcement"], mphiNeg, None, None,
                                                 eletype="Beam", oppReinf=model_pos[0]["reinforcement"], pflag=False)
        # Any warnings
        w["neg"] = [self.WARN_ELE_MAX, self.WARN_ELE_MIN]

        if model_temp is not None:
            model_neg = model_temp
            hinge_models_neg = model_temp[4]

        return model_pos, hinge_models_pos, model_neg, hinge_models_neg, w

    def model_to_df(self, model):
        """
        Main purpose of the function is to transform the hinge model dictionary into a DataFrame for use in RCMRF
        :return: DataFrame                          Lumped hinge model information for RCMRF Hysteretic model
        """
        # Initialize DataFrame with its columns
        columns = ["Element", "Bay", "Storey", "Direction", "Position", "b", "h", "coverNeg", "coverPos", "lp",
                   "length", "phi1Neg", "phi2Neg", "phi3Neg", "m1Neg", "m2Neg", "m3Neg", "phi1", "phi2", "phi3", "m1",
                   "m2", "m3"]

        numericCols = ["Bay", "Storey", "Direction", "b", "h", "coverNeg", "coverPos", "lp", "length", "phi1Neg",
                       "phi2Neg", "phi3Neg", "m1Neg", "m2Neg", "m3Neg", "phi1", "phi2", "phi3", "m1", "m2", "m3"]

        df = pd.DataFrame(columns=columns)
        for ele in model:
            if ele.lower() == "beams":
                mTemp = model[ele]["Pos"]
            else:
                mTemp = model[ele]

            for j in mTemp:
                bay = int(j[-1])
                st = int(j[-3])
                if bay == 1 or bay == self.nbays:
                    pos = "external"
                else:
                    pos = "internal"
                lp = mTemp[j][0]["lp"]

                if ele.lower() == "beams":
                    temp = np.array([ele[:-1], bay, st, self.direction, pos, model[ele]["Neg"][j][0]["b"],
                                     model[ele]["Neg"][j][0]["h"], model[ele]["Neg"][j][0]["cover"],
                                     model[ele]["Pos"][j][0]["cover"], lp,  self.bay_widths[bay-1]])
                    phiNeg = model[ele]["Neg"][j][4]["phi"][1:]
                    mNeg = model[ele]["Neg"][j][4]["m"][1:]
                    phiPos = model[ele]["Pos"][j][4]["phi"][1:]
                    mPos = model[ele]["Pos"][j][4]["m"][1:]
                else:
                    temp = np.array([ele[:-1], bay, st, self.direction, pos, model[ele][j][0]["b"],
                                     model[ele][j][0]["h"], model[ele][j][0]["cover"], model[ele][j][0]["cover"], lp,
                                     self.heights[st-1]])
                    phiNeg = phiPos = model[ele][j][4]["phi"][1:]
                    mNeg = mPos = model[ele][j][4]["m"][1:]

                data = np.concatenate((temp, phiNeg, mNeg, phiPos, mPos)).reshape(1, len(columns))

                # Concatenate into the DataFrame
                df = df.append(pd.DataFrame(data=data, columns=columns), ignore_index=True)

                # Add symmetric elements
                bayCount = self.nbays - (bay - 1) if ele.lower() == "beams" else self.nbays + 2 - bay
                if bay != bayCount:
                    bay = bayCount
                    temp[1] = bay
                    data = np.concatenate((temp, phiNeg, mNeg, phiPos, mPos)).reshape(1, len(columns))
                    # Concatenate into the DataFrame
                    df = df.append(pd.DataFrame(data=data, columns=columns), ignore_index=True)

        for i in numericCols:
            df[i] = pd.to_numeric(df[i])

        return df

    def estimate_ductilities(self, details, modes):
        """
        Estimates system hardening ductility
        :param details: dict                    Moment-curvature relationships of the elements
        :param modes: dict                      Periods and modal shapes obtained from modal analysis
        :return: float                          Hardening ductility
        """
        p = Plasticity(lp_name="Priestley", db=20, fy=self.fy, fu=self.fy * self.k_hard, lc=self.heights)
        mu_c, mu_f = p.estimate_ductilities(self.dy, details, modes)
        return mu_c, mu_f

