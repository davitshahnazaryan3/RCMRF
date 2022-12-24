import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import pickle

from design.design_utils import _create_path, get_critical_designs
from design.detailingRCMRF import DetailingRCMRF

from utils.utils import getIndex

from design.analysisMethods import run_opensees_analysis
from design.input import Input
from design.elasticAnalysis import ElasticAnalysis


class EurocodeDesign:
    def __init__(self, case, site, flag3d, soil_class="C", type_spectra=1, importance_class=2, damping=0.05,
                 fstiff=0.5):
        self.case = case
        self.site = site
        self.flag3d = flag3d
        self.output_path = _create_path(self.case)
        self.data = None

        self.soil_class = soil_class
        self.type_spectra = type_spectra
        self.importance_class = importance_class
        self.damping = damping
        self.fstiff = fstiff
        self.spectra = dict()

    def read_input(self):
        self.data = Input(self.site, self.flag3d)

    def read_hazard(self, hazard_path):
        with open(hazard_path, 'rb') as file:
            hazard = pickle.load(file)

        return hazard

    def get_preliminary_solutions(self, path):
        solution_y = None
        solution_gr = None

        # Cross-section files
        if self.flag3d:
            solution_y = pd.read_csv(path / "c-s-y.csv", index_col=0).iloc[0]
            solution_gr = pd.read_csv(path / "c-s-gr.csv", index_col=0).iloc[0]

        solution_x = pd.read_csv(path / "c-s-x.csv", index_col=0).iloc[0]
        solution = {"x_seismic": solution_x, "y_seismic": solution_y, "gravity": solution_gr}

        return solution

    def get_spectrum_parameters(self):
        soil_class = self.soil_class
        type_spectra = self.type_spectra
        damping = self.damping
        eta = max(np.sqrt(10 / (5 + damping * 100)), 0.55)
        if type_spectra == 1:
            if soil_class == 'A':
                S = 1.0
                Tb = 0.15
                Tc = 0.4
                Td = 2.
            elif soil_class == 'B':
                S = 1.2
                Tb = 0.15
                Tc = 0.5
                Td = 2.
            elif soil_class == 'C':
                S = 1.15
                Tb = 0.2
                Tc = 0.6
                Td = 2.
            elif soil_class == 'D':
                S = 1.35
                Tb = 0.2
                Tc = 0.8
                Td = 2.
            elif soil_class == 'F':
                S = 1.4
                Tb = 0.15
                Tc = 0.5
                Td = 2.
            else:
                raise ValueError("Wrong soil type selected")
        elif type_spectra == 2:
            if soil_class == 'A':
                S = 1.0
                Tb = 0.05
                Tc = 0.25
                Td = 1.2
            elif soil_class == 'B':
                S = 1.35
                Tb = 0.05
                Tc = 0.25
                Td = 1.2
            elif soil_class == 'C':
                S = 1.5
                Tb = 0.1
                Tc = 0.25
                Td = 1.2
            elif soil_class == 'D':
                S = 1.8
                Tb = 0.1
                Tc = 0.3
                Td = 1.2
            elif soil_class == 'F':
                S = 1.6
                Tb = 0.05
                Tc = 0.25
                Td = 1.2
            else:
                raise ValueError("Wrong soil type selected")
        else:
            raise ValueError("Wrong spectrum type selected")

        return S, Tb, Tc, Td, eta

    def run_modal_analysis(self, solution):
        # Create the 3D model and perform modal analysis to identify the actual period
        model_periods, modalShape, gamma, mstar = run_opensees_analysis(self.data.spans_x, self.data.spans_y,
                                                                        self.data.heights, solution,
                                                                        self.data.Ec, self.fstiff,
                                                                        self.data.inputs["seismic"])

        return model_periods, modalShape, gamma, mstar

    def get_design_gravity_loads(self):

        design_loads = [0] * self.data.nst
        seismic = [0] * self.data.nst
        for st in range(self.data.nst):
            if st == self.data.nst - 1:
                design_loads[st] = 1.35 * self.data.dead_load + 1.5 * self.data.roof_load
                seismic[st] = self.data.dead_load + 0.3 * self.data.roof_load
            else:
                design_loads[st] = 1.35 * self.data.dead_load + 1.5 * self.data.live_load
                seismic[st] = self.data.dead_load + 0.3 * self.data.live_load

        self.data.inputs["loads"] = design_loads
        self.data.inputs["seismic"] = seismic

        self.data.get_masses()

    def get_ECelastic_spectra(self, PGA):
        # elastic RS EC8 3.2.2.2
        # Type 1 spectra
        S, Tb, Tc, Td, eta = self.get_spectrum_parameters()
        T = np.linspace(0., 4., 401)
        # Sa in g, Sd in cm
        Sa = np.array([])
        Sd = np.array([])
        for i in range(len(T)):
            if T[i] <= Tb:
                Sa = np.append(Sa, (PGA * S * (1 + T[i] / Tb * (eta * 2.5 - 1))))
            elif Tb < T[i] <= Tc:
                Sa = np.append(Sa, (PGA * S * eta * 2.5))
            elif Tc < T[i] <= Td:
                Sa = np.append(Sa, (PGA * S * eta * 2.5 * Tc / T[i]))
            elif Td < T[i] <= 4:
                Sa = np.append(Sa, (PGA * S * eta * 2.5 * Tc * Td / T[i] ** 2))
            else:
                print('Wrong period range!')
            Sd = np.append(Sd, (Sa[i] * 9.81 * T[i] ** 2 / 4 / np.pi ** 2 * 100))
        return T, Sa

    def get_ECdesign_spectra(self, PGA, q, beta=0.2):
        # elastic RS EC8 3.2.2.2
        # Type 1 spectra
        S, Tb, Tc, Td, eta = self.get_spectrum_parameters()
        T = np.linspace(0., 4., 401)
        # Sa in g, Sd in cm
        Sa = np.array([])
        Sd = np.array([])
        for i in range(len(T)):
            if T[i] <= Tb:
                Sa = np.append(Sa, (PGA * S * (2 / 3 + T[i] / Tb * (2.5 / q - 2 / 3))))
            elif Tb < T[i] <= Tc:
                Sa = np.append(Sa, (PGA * S * 2.5 / q))
            elif Tc < T[i] <= Td:
                Sa = np.append(Sa, (max(beta * PGA, PGA * S * 2.5 / q * Tc / T[i])))
            elif Td < T[i] <= 4:
                Sa = np.append(Sa, (max(beta * PGA, PGA * S * 2.5 / q * Tc * Td / T[i] ** 2)))
            else:
                print('Wrong period range!')
            Sd = np.append(Sd, (Sa[i] * 9.81 * T[i] ** 2 / 4 / np.pi ** 2 * 100))
        return T, Sa

    def apply_ec_based_analysis(self, solution, model_periods, modalShape, hazard, TR=475, ductility_class="DCM"):
        # Get spectral parameters
        S, Tb, Tc, Td, eta = self.get_spectrum_parameters()

        # Get heights and number of storeys
        data = self.data
        heights = data.heights
        nst = data.nst

        # Initialize the ELFs for each direction
        if self.flag3d:
            Fi = np.zeros((2, nst))
        else:
            Fi = np.zeros((1, nst))

        hinge_models = {"x_seismic": None, "y_seismic": None, "gravity": None}
        details_models = {"x_seismic": None, "y_seismic": None, "gravity": None}
        demands_gravity = {}

        # Get masses
        masses = data.masses
        for i in range(Fi.shape[0]):
            T1 = model_periods[i]
            T1_tag = 'SA(%.2f)' % 0.3 if round(T1, 1) == 0.3 else 'SA(%.1f)' % T1

            # Computation of PGA
            Hs = hazard[2][0]
            sa_range = hazard[1][0]
            interpolation = interp1d(Hs, sa_range)
            PGA = interpolation(1 / TR)

            # Computation of SaT1
            idx = int(round(T1 * 10, 0))
            Hs = hazard[2][idx]
            sa_range = hazard[1][idx]
            interpolation = interp1d(Hs, sa_range)
            SaT1 = interpolation(1 / TR)
            print(f"[PERIOD] {T1}s, [SaT1] {SaT1}g, [PGA] {PGA}g")

            # Get elastic spectrum
            T, Sa = self.get_ECelastic_spectra(PGA)

            # EC8 table 5.1, 5.2.2.2
            q0 = 3 * 1.3 if ductility_class == 'DCM' else 4.5 * 1.3

            # for frame and frame equivalent dual systems
            kw = 1
            q = max(1.5, q0 * kw)
            if self.importance_class == 1:
                yI = 0.8
            elif self.importance_class == 2:
                yI = 1.0
            elif self.importance_class == 3:
                yI = 1.2
            else:
                yI = 1.4

            # Scale the spectrum to match the hazard
            SaFactor = float(Sa[getIndex(T1, T)] / SaT1)
            Sa = Sa / SaFactor / q if T1 >= Tb else (5 / 3 + T1 / Tb * (2.5 / q - 2 / 3)) / (
                    1 + T1 / Tb * (2.5 - 1)) * Sa / SaFactor
            self.spectra[i] = {"T": T, "Sa": Sa}

            Lam = 0.85 if (T1 <= 2 * Tc) and (nst > 2) else 1
            SaT1 = Sa[getIndex(T1, T)] * yI

            # Calculate lateral force distribution
            m = masses / data.n_seismic
            Fb = SaT1 * 9.81 * sum(m) * Lam

            print(f"[SaT1 after scaling] {SaT1:.2}, Base shear {Fb:.2}")

            z = np.cumsum(heights)
            Fi[i, :] = np.array([float(Fb * m[i] * z[i] / sum(map(lambda x, y: x * y, z, m))) for i in range(nst)])

            # Getting the demands by applying the forces along each direction sequentially
            op = ElasticAnalysis(self.data.spans_x, self.data.spans_y, self.data.heights, solution, self.data.Ec,
                                 self.fstiff, self.data.inputs, flag3d=self.flag3d)

            demands_gr, diagrams_gr = op.run_elastic_analysis(gravity=True, lateral=False, direction=i)

            # For DCL only gravity design is carried out following EC2 recommendations
            demands_elfm_gr, diagrams_elfm_gr = op.run_elastic_analysis(gravity=True, lateral=True, direction=i,
                                                                        lat_action=Fi[i, :])

            # Postprocess analysis results
            if self.flag3d:
                demands = self.postprocess_analysis_results(demands_gr, diagrams_gr, demands_elfm_gr, diagrams_elfm_gr,
                                                            direction=i)

                gr_id = "x" if i == 0 else "y"
                print(f"[DESIGN] {gr_id} direction")
                demands_gravity[gr_id] = demands["gravity"]
                seismic_demands = demands["x_seismic"] if i == 0 else demands["y_seismic"]
                seismic_solution = solution["x_seismic"] if i == 0 else solution["y_seismic"]

            else:
                demands = demands_elfm_gr
                seismic_demands = demands
                seismic_solution = solution["x_seismic"]

            modes = modalShape[:, 0]

            # Detail the structural elements
            details, hinges, mu_c, mu_f, warnMax, warnMin, warnings = \
                self.design_elements(seismic_demands, seismic_solution, modes, None, cover=0.03,
                                     direction=i, est_ductilities=False)

            warnings_min = sum(warnings["MIN"]["Columns"].values())
            warnings_max = sum(warnings["MAX"]["Columns"].values())
            n_elements = len(warnings["MIN"]["Columns"])

            print(f"[WARNINGS] MIN: {warnings_min / n_elements * 100}%")
            print(f"[WARNINGS] MAX: {warnings_max / n_elements * 100}%")

            if i == 0:
                hinge_models["x_seismic"] = hinges
                details_models["x_seismic"] = details
            else:
                hinge_models["y_seismic"] = hinges
                details_models["y_seismic"] = details

        if self.flag3d:
            # Design the interior elements (envelope of both directions)
            details, hinge_gravity, warn_gravity = self.design_elements(demands_gravity, solution["gravity"], None, None,
                                                                        cover=0.03, direction=0, gravity=True,
                                                                        est_ductilities=False)
            warnings_min = sum(warn_gravity["warnings"]["MIN"]["Columns"].values())
            warnings_max = sum(warn_gravity["warnings"]["MAX"]["Columns"].values())
            n_elements = len(warn_gravity["warnings"]["MIN"]["Columns"])
            print(f"[WARNINGS] MIN: {warnings_min / n_elements * 100}%")
            print(f"[WARNINGS] MAX: {warnings_max / n_elements * 100}%")
            hinge_models["gravity"] = hinge_gravity
            details_models["gravity"] = details

            hinge_x, hinge_y = get_critical_designs(hinge_models["x_seismic"], hinge_models["y_seismic"])
            hinge_models["x_seismic"] = hinge_x
            hinge_models["y_seismic"] = hinge_y

        return hinge_models, details_models

    def design_elements(self, demands, sections, modes, dy, ductility_class="DCM", cover=0.03, est_ductilities=True,
                        direction=0, gravity=False):
        """
        Runs M-phi to optimize for reinforcement for each section
        :param demands: DataFrame or dict           Demands identified from a structural analysis (lateral+gravity)
        :param sections: DataFrame                  Solution including section information
        :param modes: dict                          Periods and modal shapes obtained from modal analysis
        :param dy: float                            System yield displacement in m
        :param ductility_class: str                 Ductility class (DCM or DCH, following Eurocode 8 recommendations)
        :param cover: float                         Reinforcement cover in m
        :param est_ductilities: bool                Whether to estimate hardening and fracturing ductilities
        :param direction: bool                      0 for x direction, 1 for y direction
        :param gravity: bool                        Design gravity frames condition
        :return: dict                               Designed element properties from the moment-curvature relationship
        """
        if direction == 0:
            nbays = self.data.n_bays
            spans = self.data.spans_x
        else:
            spans = self.data.spans_y
            nbays = len(spans)

        d = DetailingRCMRF(demands, self.data.nst, nbays, self.data.fy, self.data.fc, spans, self.data.heights,
                           self.data.n_seismic, self.data.masses, dy, sections, ductility_class=ductility_class,
                           rebar_cover=cover, est_ductilities=est_ductilities, direction=direction)
        if gravity:
            data, hinge_models, w = d.design_gravity()
            warnMax = d.WARNING_MAX
            warnMin = d.WARNING_MIN
            warnings = {"warnings": w, "warnMax": warnMax, "warnMin": warnMin}

            return data, hinge_models, warnings

        else:
            data, hinge_models, mu_c, mu_f, warnings = d.design_elements(modes)
            warnMax = d.WARNING_MAX
            warnMin = d.WARNING_MIN

        return data, hinge_models, mu_c, mu_f, warnMax, warnMin, warnings

    def postprocess_analysis_results(self, gravity, diagrams_gravity, elfm=None, diagrams_elfm=None, direction=0):
        if direction == 0:
            # Along X axis
            midx = 5
        else:
            # Along Y axis
            midx = 4

        if elfm:
            # Get critical of ELFM and ELFM+gravity
            for key in gravity:
                for ele in gravity[key]:
                    for force in gravity[key][ele]:
                        if ele.startswith("Beams") and force == "M":
                            # Beams and M
                            for d in gravity[key][ele][force]:
                                """
                                Note: fair choice is selection of corresponding values of other internal forces, rather
                                than the max for each, but for now, we stick with the max values for each force
                                """
                                gravity[key][ele][force][d] = np.maximum(gravity[key][ele][force][d],
                                                                         elfm[key][ele][force][d])

                        else:
                            gravity[key][ele][force] = np.maximum(gravity[key][ele][force],
                                                                  elfm[key][ele][force])

        # Get the max internal forces from diagrams of beams
        # Diagrams - 0-N, 1-Vy, 2-Vz, 3-T, 4-My, 5-Mz
        for ele in diagrams_gravity:
            if diagrams_elfm:
                critical_pos = min(min(diagrams_gravity[ele][:, midx]), min(diagrams_elfm[ele][:, midx]))
                critical_neg = max(max(diagrams_gravity[ele][:, midx]), max(diagrams_elfm[ele][:, midx]))
            else:
                critical_pos = min(diagrams_gravity[ele][:, midx])
                critical_neg = max(diagrams_gravity[ele][:, midx])

            ele = str(ele)
            if ele[0] == "3":
                if 1 < int(ele[2]) < len(self.data.spans_y):
                    # Gravity, x direction
                    bay = int(ele[1]) - 1
                    st = int(ele[3]) - 1
                    ybay = int(ele[2]) - 1

                    gravity = self.assign_demands(gravity, "gravity", critical_pos, critical_neg, bay, st,
                                                  add_tag="Beams_x", ybay=ybay)

                else:
                    # x_seismic
                    bay = int(ele[1]) - 1
                    st = int(ele[3]) - 1
                    gravity = self.assign_demands(gravity, "x_seismic", critical_pos, critical_neg, bay, st)

            else:
                if 1 < int(ele[1]) < self.data.n_bays:
                    # Gravity, y direction
                    bay = int(ele[1]) - 1
                    st = int(ele[3]) - 1
                    ybay = int(ele[2]) - 1

                    gravity = self.assign_demands(gravity, "gravity", critical_pos, critical_neg, bay, st,
                                                  add_tag="Beams_y", ybay=ybay)

                else:
                    # y_seismic
                    bay = int(ele[2]) - 1
                    st = int(ele[3]) - 1
                    gravity = self.assign_demands(gravity, "y_seismic", critical_pos, critical_neg, bay, st)

        return gravity

    @staticmethod
    def assign_demands(demands, tag, pos, neg, bay, st, add_tag="Beams", ybay=None):
        if tag == "gravity":
            if pos > demands[tag][add_tag]["M"]["Pos"][st][bay][ybay]:
                demands[tag][add_tag]["M"]["Pos"][st][bay][ybay] = pos
            if neg > demands[tag][add_tag]["M"]["Neg"][st][bay][ybay]:
                demands[tag][add_tag]["M"]["Pos"][st][bay][ybay] = neg
        else:
            if pos > demands[tag][add_tag]["M"]["Pos"][st][bay]:
                demands[tag][add_tag]["M"]["Pos"][st][bay] = pos
            if neg > demands[tag][add_tag]["M"]["Neg"][st][bay]:
                demands[tag][add_tag]["M"]["Pos"][st][bay] = neg

        return demands
