"""
Defines lumped hinge models
NOTE: Important assumption - except for Haselton springs, everywhere else (e.g. modelling the elements, Ec is taken
equal to the value submitted from the IPBSD framework (i.e. cracked concrete properties).
For Haselton springs, gross EI is utilized, as the subsequent expressions for effective stiffness values account for all
contributors to pre-yield displacements, including flexure and shear strains, concrete cracking, and bond-slip.
Ref: Haselton et al. (2016) Calibration of model to simulate response of reinforced concrete beam-columns to collapse.
"""
import openseespy.opensees as op


class Sections:

    def __init__(self, sections, materials):
        """
        Initializes hinge model creation
        :param sections:  DataFrame                 Element properties
        :param materials: dict                      Concrete and reinforcement material properties
        """
        self.sections = sections
        self.materials = materials
        self.UBIG = 1.0e10

    def _rot_spring_2d_modIKmodel(self, eleID, nodeR, nodeC, K, asPos, asNeg, MyPos, MyNeg, LS, LK, LA, LD, cS, cK, cA,
                                  cD, th_pP, th_pN, th_pcP, th_pcN, ResP, ResN, th_uP, th_uN, DP, DN):
        """
        This routine creates a uniaxial material spring with deterioration
        Spring follows: Bilinear Response based on Modified Ibarra Krawinkler Deterioration Model
        Written by: Dimitrios G. Lignos, Ph.D.
        :param eleID: int                           Element identification
        :param nodeR: int                           Retained/master node
        :param nodeC: int                           Constrained/slave node
        :param K: float                             Initial stiffness after the modification for n (see Ibarra and
                                                                                                    Krawinkler, 2005)
        :param asPos: float                         Strain hardening ratio after n modification (see Ibarra and
                                                                                                    Krawinkler, 2005)
        :param asNeg: float                         Strain hardening ratio after n modification (see Ibarra and
                                                                                                    Krawinkler, 2005)
        :param MyPos: float                         Positive yield moment (with sign)
        :param MyNeg: float                         Negative yield moment (with sign)
        :param LS: float                            Basic strength deterioration parameter (see Lignos and
                                                                                                    Krawinkler, 2009)
        :param LK: float                            Unloading stiffness deterioration parameter (see Lignos and
                                                                                                    Krawinkler, 2009)
        :param LA: float                            Accelerated reloading stiffness deterioration parameter
                                                                                    (see Lignos and Krawinkler, 2009)
        :param LD: float                            Post-capping strength deterioration parameter
                                                                                    (see Lignos and Krawinkler, 2009)
        :param cS: float                            Exponent for basic strength deterioration
        :param cK: float                            Exponent for unloading stiffness deterioration
        :param cA: float                            Exponent for accelerated reloading stiffness deterioration
        :param cD: float                            Exponent for post-capping strength deterioration
        :param th_pP: float                         Plastic rotation capacity for positive loading direction
        :param th_pN: float                         Plastic rotation capacity for negative loading direction
        :param th_pcP: float                        Post-capping rotation capacity for positive loading direction
        :param th_pcN: float                        Post-capping rotation capacity for negative loading direction
        :param ResP: float                          Residual strength ratio for positive loading direction
        :param ResN: float                          Residual strength ratio for negative loading direction
        :param th_uP: float                         Ultimate rotation capacity for positive loading direction
        :param th_uN: float                         Ultimate rotation capacity for negative loading direction
        :param DP: float                            Rate of cyclic deterioration for positive loading direction
        :param DN: float                            Rate of cyclic deterioration for negative loading direction
        :return: None

        References:
            Ibarra, L. F., and Krawinkler, H. (2005). �Global collapse of frame structures under seismic excitations,� Technical Report 152, The John A. Blume Earthquake Engineering Research Center, Department of Civil Engineering, Stanford University, Stanford, CA.
            Ibarra, L. F., Medina, R. A., and Krawinkler, H. (2005). �Hysteretic models that incorporate strength and stiffness deterioration,� International Journal for Earthquake Engineering and Structural Dynamics, Vol. 34, No.12, pp. 1489-1511.
            Lignos, D. G., and Krawinkler, H. (2010). �Deterioration Modeling of Steel Beams and Columns in Support to Collapse Prediction of Steel Moment Frames�, ASCE, Journal of Structural Engineering (under review).
            Lignos, D. G., and Krawinkler, H. (2009). �Sidesway Collapse of Deteriorating Structural Systems under Seismic Excitations,� Technical Report 172, The John A. Blume Earthquake Engineering Research Center, Department of Civil Engineering, Stanford University, Stanford, CA.
        """
        op.uniaxialMaterial('Bilin', eleID, K, asPos, asNeg, MyPos, MyNeg, LS, LK, LA, LD, cS, cK, cA, cD, th_pP, th_pN,
                            th_pcP, th_pcN, ResP, ResN, th_uP, th_uN, DP, DN)
        op.element('zeroLength', eleID, nodeR, nodeC, '-mat', eleID, '-dir', 6)
        op.equalDOF(nodeR, nodeC, 1, 2)

    def haselton_springs(self, idx, tag=None, nodeR=None, base_spring=False):
        """
        Calculates parameters for Haselton springs
        :param idx: int                             ID of the element of interest
        :param tag: int                             Tag for creating subsequent node hinges
        :param nodeR: int                           Retained/master node ID
        :param base_spring: bool                    Calculate for base springs or not
        :return: None
        References:
            Haselton CB, Deierlein G.G. Assessing seismic collapse safety of modern reinforced concrete moment frame
            buildings. Stanford University, 2007.
        """
        MyPos = self.sections['MyPos'][idx]
        MyNeg = -self.sections['MyNeg'][idx]
        P = self.sections['Ptotal'][idx]
        Ec = (3320 * float(self.materials['fc']) ** 0.5 + 6900) * 1000.0
        b = self.sections['b'][idx]
        h = self.sections['h'][idx]
        cover_pos = self.sections['coverPos'][idx]
        cover_neg = self.sections['coverNeg'][idx]
        d_pos = h - cover_pos
        d_neg = h - cover_neg
        length = self.sections['Length'][idx]
        Ig = b * h ** 3 / 12
        Ag = b * h
        fc_prime = float(self.materials['fc']) * 1000.0
        EIg = Ec * Ig
        nu = P / Ag / fc_prime
        EIy = 0.75 * (0.1 + nu) ** 0.8 * EIg
        if EIy / EIg < 0.2:
            EIy = 0.2 * EIg
        elif EIy / EIg > 0.6:
            EIy = 0.6 * EIg
        else:
            EIy = 0.75 * (0.1 + nu) ** 0.8 * EIg
        asl = self.sections['asl'][idx]
        Ash = self.sections['Ash'][idx]
        spacing = self.sections['spacing'][idx]
        rosh = Ash / spacing / b
        db = self.sections['db'][idx] / 1000
        sn = spacing / db
        ro_pos = self.sections['ro_long_pos'][idx]
        ro_neg = self.sections['ro_long_neg'][idx]
        theta_cap_plPos = 0.12 * (1 + 0.55 * asl) * 0.16 ** nu * (0.02 + 40 * rosh) ** 0.43 * 0.54 ** \
                          (0.01 * fc_prime / 1000) * 0.66 ** (0.1 * sn) * 2.27 ** (10.0 * ro_pos)
        theta_cap_totPos = 0.12 * (1 + 0.4 * asl) * 0.2 ** nu * (0.02 + 40 * rosh) ** 0.52 * 0.56 ** \
                           (0.01 * fc_prime / 1000) * 2.37 ** (10.0 * ro_pos)
        theta_cap_plNeg = 0.12 * (1 + 0.55 * asl) * 0.16 ** nu * (0.02 + 40 * rosh) ** 0.43 * 0.54 ** \
                          (0.01 * fc_prime / 1000) * 0.66 ** (0.1 * sn) * 2.27 ** (10.0 * ro_neg)
        theta_cap_totNeg = 0.12 * (1 + 0.4 * asl) * 0.2 ** nu * (0.02 + 40 * rosh) ** 0.52 * 0.56 ** \
                           (0.01 * fc_prime / 1000) * 2.37 ** (10.0 * ro_neg)
        theta_pc = min(0.76 * 0.031 ** nu * (0.02 + 40 * rosh) ** 1.02, 0.1)
        mc_my = 1.13
        c = self.sections['c'][idx]
        D = self.sections['D'][idx]
        gamma = 30 * 0.3 ** nu
        k0 = 6 * EIy / length * 11
        res_strength = self.sections['Res'][idx]
        as_pos = 11 * (MyPos * (mc_my - 1)) / (k0 * theta_cap_plPos) / (1 + 10 * (1 - 11 * (MyPos * (mc_my - 1)) /
                                                                                  (k0 * theta_cap_plPos)))
        as_neg = 11 * (-MyNeg * (mc_my - 1)) / (k0 * theta_cap_plNeg) / (1 + 10 * (1 - 11 * (-MyNeg * (mc_my - 1)) /
                                                                                   (k0 * theta_cap_plNeg)))
        if base_spring:
            eleID = 10001 + idx
            try:
                nodeC = int(f"{nodeR}0")
                self._rot_spring_2d_modIKmodel(eleID, nodeR, nodeC, k0, as_pos, as_neg, MyPos, MyNeg, gamma, gamma, 0.0,
                                               0.0, c, c, c, c, theta_cap_plPos, theta_cap_plNeg, theta_pc, theta_pc,
                                               res_strength, res_strength, theta_cap_totPos + theta_pc,
                                               theta_cap_totNeg + theta_pc, D, D)
            except TypeError:
                print('[EXCEPTION] Master node ID not provided')

        else:
            try:
                eleID = 100000 + tag
                op.uniaxialMaterial('Bilin', eleID, k0, as_pos, as_neg, MyPos, MyNeg, gamma, gamma, gamma, gamma, c, c,
                                    c, c, theta_cap_plPos, theta_cap_plNeg, theta_pc, theta_pc, res_strength,
                                    res_strength, theta_cap_totPos + theta_pc, theta_cap_totNeg + theta_pc, D, D)
            except TypeError:
                print('[EXCEPTION] Node ID not provided')

    def hysteretic_hinges(self, et, iNode, jNode, ele, transfTag, flag3d, tcl_file):
        """
        Creates hysteretic hinges
        :param et: int                              Element tag
        :param ele: DataFrame                       Hinge model parameters
        :param transfTag: int                       Element transformation tag
        :param flag3d: bool                         True for 3D modelling, False for 2D modelling
        :return: None
        """
        # Cross-section area of the element
        area = ele['b'] * ele['h']
        # Moment of inertia of the cross-section
        iz = ele['b'] * ele['h'] ** 3 / 12
        # Secondary moment of inertia
        iy = ele['h'] * ele['b'] ** 3 / 12
        # Shear parameters
        nu = 0.2
        Gc = float(self.materials['Ec']) * 1000.0 / 2.0 / (1 + nu)
        # Torsional moment of inertia
        if ele["h"] >= ele["b"]:
            J = ele["b"] * ele["h"]**3 * (16 / 3 - 3.36 * ele["h"] / ele["b"] * (1 - 1 / 12 * (ele["h"] / ele["b"])**4))
        else:
            J = ele["h"] * ele["b"]**3 * (16 / 3 - 3.36 * ele["b"] / ele["h"] * (1 - 1 / 12 * (ele["b"] / ele["h"])**4))

        # Node IDs connecting the elements
        if not flag3d:
            # Bay and storey levels
            bay = ele['Bay']
            st = ele['Storey']

            if ele['Element'].lower() == 'beam':
                iNode = int(f"{bay}{st}")
                jNode = int(f"{bay + 1}{st}")
            else:
                iNode = int(f"{bay}{st - 1}")
                jNode = int(f"{bay}{st}")

        # Material tags
        matTag1 = int(f'101{et}')
        matTag2 = int(f'102{et}')
        intTag = int(f'105{et}')
        phTag1 = int(f'106{et}')
        phTag2 = int(f'107{et}')

        # Integration tag
        integrationTag = int(f'108{et}')

        # Some additional parameters for the hysteretic model
        pinchX = 0.8
        pinchY = 0.5
        damage1 = 0.0
        damage2 = 0.0
        beta = 0.0

        # Creating the hinges at both ends (equivalent assumption)
        op.uniaxialMaterial('Hysteretic', matTag1, ele['m1'], ele['phi1'], ele['m2'], ele['phi2'], ele['m3'],
                            ele['phi3'], -ele['m1Neg'], -ele['phi1Neg'], -ele['m2Neg'], -ele['phi2Neg'], -ele['m3Neg'],
                            -ele['phi3Neg'], pinchX, pinchY, damage1, damage2, beta)
        op.uniaxialMaterial('Hysteretic', matTag2, ele['m1'], ele['phi1'], ele['m2'], ele['phi2'], ele['m3'],
                            ele['phi3'], -ele['m1Neg'], -ele['phi1Neg'], -ele['m2Neg'], -ele['phi2Neg'], -ele['m3Neg'],
                            -ele['phi3Neg'], pinchX, pinchY, damage1, damage2, beta)

        # Writing to tcl
        if tcl_file:
            tcl_file.write(f"\nuniaxialMaterial Hysteretic {matTag1} {ele['m1']} {ele['phi1']} {ele['m2']} "
                           f"{ele['phi2']} {ele['m3']} {ele['phi3']} -{ele['m1Neg']} -{ele['phi1Neg']} -{ele['m2Neg']}"
                           f" -{ele['phi2Neg']} -{ele['m3Neg']} {-ele['phi3Neg']} {pinchX} {pinchY} {damage1}"
                           f" {damage2} {beta};")
            tcl_file.write(f"\nuniaxialMaterial Hysteretic {matTag2} {ele['m1']} {ele['phi1']} {ele['m2']} "
                           f"{ele['phi2']} {ele['m3']} {ele['phi3']} -{ele['m1Neg']} -{ele['phi1Neg']} -{ele['m2Neg']}"
                           f" -{ele['phi2Neg']} -{ele['m3Neg']} {-ele['phi3Neg']} {pinchX} {pinchY} {damage1}"
                           f" {damage2} {beta};")

        # Elastic section
        if flag3d:
            op.section('Elastic', intTag, float(self.materials['Ec']) * 1000.0, area, iy, iz, Gc, J)
            if tcl_file:
                tcl_file.write(f"\nsection Elastic {intTag} {float(self.materials['Ec']) * 1000.0} {area} {iy} {iz} "
                               f"{Gc} {J};")
        else:
            op.section('Elastic', intTag, float(self.materials['Ec']) * 1000.0, area, iz)

        # Create the plastic hinge flexural section about ZZ
        op.section('Uniaxial', phTag1, matTag1, 'Mz')
        op.section('Uniaxial', phTag2, matTag2, 'Mz')

        if tcl_file:
            tcl_file.write(f"\nsection Uniaxial {phTag1} {matTag1} Mz;")
            tcl_file.write(f"\nsection Uniaxial {phTag2} {matTag2} Mz;")

        aggTag1 = int(f"114{et}")
        aggTag2 = int(f"115{et}")
        if transfTag == 1 and flag3d:
            # Transformation tag 1 refers to the columns
            # Additional hinges are required for bidirectional response for the columns in the 3D model
            # Since those are symmetrical square columns, i.e. designed to have the same properties along both principal
            # directions, then the properties will match.
            # Additional tags for the materials
            matTag3 = int(f'111{et}')
            matTag4 = int(f'112{et}')
            axialTag = int(f'113{et}')

            # Beam integration
            op.beamIntegration('HingeRadau', integrationTag, aggTag1, ele['lp'], aggTag2, ele['lp'], intTag)

            # Create the plastic hinge axial material
            op.uniaxialMaterial("Elastic", axialTag, float(self.materials['Ec']) * 1000.0 * area)

            # Create the plastic hinge materials
            op.uniaxialMaterial('Hysteretic', matTag3, ele['m1'], ele['phi1'], ele['m2'], ele['phi2'], ele['m3'],
                                ele['phi3'], -ele['m1Neg'], -ele['phi1Neg'], -ele['m2Neg'], -ele['phi2Neg'],
                                -ele['m3Neg'], -ele['phi3Neg'], pinchX, pinchY, damage1, damage2, beta)
            op.uniaxialMaterial('Hysteretic', matTag4, ele['m1'], ele['phi1'], ele['m2'], ele['phi2'], ele['m3'],
                                ele['phi3'], -ele['m1Neg'], -ele['phi1Neg'], -ele['m2Neg'], -ele['phi2Neg'],
                                -ele['m3Neg'], -ele['phi3Neg'], pinchX, pinchY, damage1, damage2, beta)

            # Aggregate P and Myy behaviour to Mzz behaviour
            op.section("Aggregator", aggTag1, axialTag, "P", matTag3, "My", "-section", phTag1)
            op.section("Aggregator", aggTag2, axialTag, "P", matTag4, "My", "-section", phTag2)

            if tcl_file:
                tcl_file.write(f"\nuniaxialMaterial Elastic {axialTag} {float(self.materials['Ec']) * 1000.0 * area};")
                tcl_file.write(f"\nuniaxialMaterial Hysteretic {matTag3} {ele['m1']} {ele['phi1']} {ele['m2']} "
                               f"{ele['phi2']} {ele['m3']} {ele['phi3']} -{ele['m1Neg']} -{ele['phi1Neg']}"
                               f" -{ele['m2Neg']} -{ele['phi2Neg']} -{ele['m3Neg']} {-ele['phi3Neg']} "
                               f"{pinchX} {pinchY} {damage1} {damage2} {beta};")
                tcl_file.write(f"\nuniaxialMaterial Hysteretic {matTag4} {ele['m1']} {ele['phi1']} {ele['m2']} "
                               f"{ele['phi2']} {ele['m3']} {ele['phi3']} -{ele['m1Neg']} -{ele['phi1Neg']}"
                               f" -{ele['m2Neg']} -{ele['phi2Neg']} -{ele['m3Neg']} {-ele['phi3Neg']} "
                               f"{pinchX} {pinchY} {damage1} {damage2} {beta};")

                tcl_file.write(f"\nsection Aggregator {aggTag1} {axialTag} P {matTag3} My -section {phTag1};")
                tcl_file.write(f"\nsection Aggregator {aggTag2} {axialTag} P {matTag4} My -section {phTag2};")

                lp = round(ele['lp'], 3)
                tcl_file.write(f'\nelement forceBeamColumn {int(et)} {iNode} {jNode} {transfTag} "HingeRadau {aggTag1}'
                               f' {lp} {aggTag2} {lp} {intTag}";')

        else:
            # Beam integration
            op.beamIntegration('HingeRadau', integrationTag, phTag1, ele['lp'], phTag2, ele['lp'], intTag)

            lp = round(ele['lp'], 3)

            if tcl_file:
                tcl_file.write(f'\nelement forceBeamColumn {int(et)} {iNode} {jNode} {transfTag} "HingeRadau {phTag1}'
                               f' {lp} {phTag2} {lp} {intTag}";')

        op.element('forceBeamColumn', int(et), iNode, jNode, transfTag, integrationTag)
