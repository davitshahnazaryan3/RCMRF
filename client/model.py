"""
Model creator of an RC MRF. Lumped hinge models following the recommendations of Haselton 2007 are used.
"""
import os
import openseespy.opensees as op
import pandas as pd, json, pickle
from client.geometry import Geometry
from client.sections import Sections
from client.recorders import Recorders
from analysis.static import Static
from analysis.modal import Modal
from analysis.spo import SPO


class Model:

    def __init__(self, analysis_type, sections_file, loads_file, materials, outputsDir, system='Perimeter',
                 hingeModel='haselton'):
        """
        Initializes OpenSees model creator
        :param analysis_type: list(str)             Type of analysis for which we are recording [TH, PO, ST, MA, ELF]
                                                    TH - time history analysis
                                                    PO - static pushover analysis
                                                    ST - static analysis (e.g. gravity)
                                                    MA - modal analysis
                                                    ELF - equivalent lateral force method of analysis
        :param sections_file: str                   Name of the input csv file containing information on the elements
            Features:   Element - element type (column/beam)
                        Bay - bay location of element counting from the left side
                        Storey - storey location of element
                        Position - external or internal
                        b - element width in m
                        h - element height in m
                        cover - reinforcement cover in m
                        Pgrav - axial force from gravity loads in kN
                        Plateral - axial force from lateral loads in kN
                        MyPos - positive bending moment capacity in kNm
                        MyNeg - negative bending moment capacity in kNm
                        asl - bond-slip indicator variable (1=slip possible, 0=slip is not possible)
                        Ash - area of transverse reinforcement in m2
                        spacing - spacing of transverse reinforcement in m
                        db - longitudinal reinforcement diameter in mm
                        c - deterioration exponents
                        D - rate of cyclic deterioration
                        Res - residual strength
                        Length - length of element in m
                        ro_long - longitudinal reinforcement ratio
        :param loads_file: str                      Name of the csv file containing information on the loads and masses
            Features:   Storey - storey level of load
                        Pattern - distributed or point internal or point external or pdelta (necessary for perimeter
                                frames) or mass
                        Load - load value in kN/m in positive sign for 'distributed' loads or kN for 'point internal'
                            or 'point external' or 'pdelta' loads or ton for 'mass'
        :param materials: dict                      Material properties of reinforcement and concrete
        :param system: str                          MRF type, i.e. Perimeter or space
        :param hingeModel: str                      Hinge model type (Haselton (4 nodes) or Hysteretic (2 nodes))
        """
        self.base_nodes = None
        self.base_cols = None
        self.elements = None
        self.analysis_type = analysis_type
        self.materials = pd.read_csv(materials)
        self.system = system
        self.hingeModel = hingeModel.lower()
        self.COL_TRANSF_TAG = 1
        self.BEAM_TRANSF_TAG = 2
        self.NEGLIGIBLE = 1e-09
        self.outputsDir = outputsDir
        self.sections = pd.read_csv(sections_file)
        cols = [i for i in self.sections.columns if i not in ('Element', 'Position',
                                                              'Storey', 'Bay')]
        for col in cols:
            self.sections[col] = self.sections[col].astype(float)
        else:
            self.loads = pd.read_csv(loads_file)
            self.check_integrity()
            self.g = Geometry(self.sections, self.hingeModel)
            self.NUM_MODES = 3
            self.DAMP_MODES = [1, 3]
            self.results = {}

    def createFolder(self, directory):
        """
        Checks whether provided directory exists, if no creates one
        :param directory: str
        :return: None
        """
        try:
            if not os.path.exists(directory):
                os.makedirs(directory)
        except OSError:
            print('Error: Creating directory. ' + directory)

    def check_integrity(self):
        """
        Checks whether the input arguments have been supplied properly
        :return: None
        """
        if self.system not in ('Perimeter', 'Space'):
            raise ValueError('[EXCEPTION] Wrong system type provided, should be Perimeter or Space')
        print('[SUCCESS] Integrity of input arguments checked')

    def create_model(self):
        """
        Initiates model creation
        :return: None
        """
        op.model('Basic', '-ndm', 2, '-ndf', 3)
        print('[INITIATE] Model generation started')

    def create_nodes(self, fixity='fixed'):
        """
        Creates nodes
        :param fixity: str                          Boundary condition of base nodes
        :return: list                               Base node IDs
        """
        if fixity == 'fixed':
            fix = 1
        elif fixity == 'pinned':
            fix = 0
        else:
            raise Exception('[EXCEPTION] Wrong boundary condition for base nodes')

        df = self.g.define_nodes()
        base_nodes = []
        hinge_nodes = []
        for n in df.index:
            op.node(int(df['Node id'][n]), df['x'][n], df['z'][n])
            if df['z'][n] == 0 and int(df['Node id'][n]) < 10000:
                op.fix(int(df['Node id'][n]), 1, 1, fix)
                base_nodes.append(int(df['Node id'][n]))
            elif 0.0 <= df['x'][n] <= max(self.g.widths):
                if int(df['Node id'][n]) < 10000 and df['z'][n] <= max(self.g.heights):
                    hinge_nodes.append(int(df['Node id'][n]))

        print('[SUCCESS] Geometric properties have been defined')

        return base_nodes, hinge_nodes

    def define_transformations(self, col_transf_type='PDelta', beam_transf_tag='Linear'):
        """
        Defines geometric transformations for beams and columns (PDelta, Linear, Corotational)
        :param col_transf_type: str                 Column transformation type
        :param beam_transf_tag: str                 Beam transformation type
        :return: None
        """
        if any((tag not in ('PDelta', 'Linear', 'Corotational') for tag in (col_transf_type, beam_transf_tag))):
            raise Exception('[EXCEPTION] Wrong transformation type provided')
        else:
            op.geomTransf(col_transf_type, self.COL_TRANSF_TAG)
            op.geomTransf(beam_transf_tag, self.BEAM_TRANSF_TAG)
        print('[SUCCESS] Material Properties have been defined')

    def joint_materials(self):
        """
        Defines joint materials
        :return: None
        """
        num_hinge = (len(self.g.heights) - 1) * len(self.g.widths)
        Ec = float(self.materials['Ec'])
        for h in range(num_hinge):
            op.uniaxialMaterial('Elastic', 200000 + h, Ec * 1000)

        print('[SUCCESS] Joint material properties have been defined')

    def rot_springs(self, base_nodes):
        """
        Defines rotational springs at the base of the structure
        :param base_nodes: list                     Base node IDs
        :return: None
        """
        sections = self.sections[((self.sections['Element'] == 'Column') &
                                  (self.sections['Storey'] == 1))].reset_index(drop=True)
        s = Sections(sections, self.materials)
        for base_node in range(self.g.nbays + 1):
            s.haselton_springs(base_node, nodeR=(base_nodes[base_node]), base_spring=True)

        print('[SUCCESS] Rotational springs have been defined')

    def bilin_springs(self, nodes):
        """
        Defines bilinear springs for the plastic hinges of the structure
        :param nodes: list                          List of the node IDs
        :return: None
        """

        def get_last_digits(num, check_list, last_digits_count=2):
            if int(str(num)[-last_digits_count:]) in check_list:
                return num

        sections_beams = self.sections[(self.sections['Element'] == 'Beam')].reset_index(drop=True)
        sections_cols = self.sections[(self.sections['Element'] == 'Column')].reset_index(drop=True)
        beam_nodes = list(filter(None, map(lambda x: get_last_digits(x, check_list=[20, 40]), nodes)))
        col_nodes = list(filter(None, map(lambda x: get_last_digits(x, check_list=[10, 30]), nodes)))
        for n in range(len(beam_nodes)):
            storey = int(str(beam_nodes[n])[0]) - 1
            bay = int(str(beam_nodes[n])[1]) if str(beam_nodes[n])[-2:] == '20' else int(str(beam_nodes[n])[1]) - 1
            sections = sections_beams[((sections_beams['Bay'] == bay) &
                                       (sections_beams['Storey'] == storey))].reset_index(drop=True)
            s = Sections(sections, self.materials)
            s.haselton_springs(0, tag=(int(str(beam_nodes[n]))))
        else:
            for n in range(len(col_nodes)):
                storey = int(str(col_nodes[n])[0]) - 1 if str(col_nodes[n])[-2:] == '10' else int(str(col_nodes[n])[0])
                bay = int(str(col_nodes[n])[1])
                sections = sections_cols[((sections_cols['Bay'] == bay) &
                                          (sections_cols['Storey'] == storey))].reset_index(drop=True)
                s = Sections(sections, self.materials)
                s.haselton_springs(0, tag=(int(str(col_nodes[n]))))

        print('[SUCCESS] Bilinear springs have been defined')

    def create_elements(self):
        """
        Creates elastic beam column elements
        :return: dict                                   Dictionary containing all element IDs
        """
        young_modulus = float(self.materials['Ec']) * 1000.0
        elements = {'Columns internal': {},  'Columns external': {},  'Beams': {}}
        base_cols = []
        for st in range(1, self.g.nst + 1):
            elements['Columns internal'][st] = []
            elements['Columns external'][st] = []
            elements['Beams'][st] = []

        for ele in range(len(self.sections)):
            if (self.sections['Storey'][ele] == 1) & (self.sections['Element'][ele] == 'Column'):
                eleid = int(f"1{self.sections['Storey'][ele]}{self.sections['Bay'][ele]}")
                node_i = int(f"1{self.sections['Bay'][ele]}000")
                node_j = int(f"{self.sections['Storey'][ele] + 1}{self.sections['Bay'][ele]}10")
                area = self.sections['b'][ele] * self.sections['h'][ele]
                inertia = self.sections['b'][ele] * self.sections['h'][ele] ** 3 / 12
                op.element('elasticBeamColumn', eleid, node_i, node_j, area, young_modulus, inertia, self.COL_TRANSF_TAG)
                if self.sections['Bay'][ele] in (1, self.g.nbays + 1):
                    elements['Columns external'][self.sections['Storey'][ele]].append(eleid)
                else:
                    elements['Columns internal'][self.sections['Storey'][ele]].append(eleid)
                base_cols.append(eleid)
            elif (self.sections['Storey'][ele] > 1) & (self.sections['Element'][ele] == 'Column'):
                eleid = int(f"1{self.sections['Storey'][ele]}{self.sections['Bay'][ele]}")
                node_i = int(f"{self.sections['Storey'][ele]}{self.sections['Bay'][ele]}30")
                node_j = int(f"{self.sections['Storey'][ele] + 1}{self.sections['Bay'][ele]}10")
                area = self.sections['b'][ele] * self.sections['h'][ele]
                inertia = self.sections['b'][ele] * self.sections['h'][ele] ** 3 / 12
                op.element('elasticBeamColumn', eleid, node_i, node_j, area, young_modulus, inertia, self.COL_TRANSF_TAG)
                if self.sections['Bay'][ele] in (1, self.g.nbays + 1):
                    elements['Columns external'][self.sections['Storey'][ele]].append(eleid)
                else:
                    elements['Columns internal'][self.sections['Storey'][ele]].append(eleid)
            else:
                eleid = int(f"2{self.sections['Storey'][ele]}{self.sections['Bay'][ele]}")
                node_i = int(f"{self.sections['Storey'][ele] + 1}{self.sections['Bay'][ele]}20")
                node_j = int(f"{self.sections['Storey'][ele] + 1}{self.sections['Bay'][ele] + 1}40")
                area = self.sections['b'][ele] * self.sections['h'][ele]
                inertia = self.sections['b'][ele] * self.sections['h'][ele] ** 3 / 12
                op.element('elasticBeamColumn', eleid, node_i, node_j, area, young_modulus, inertia, self.BEAM_TRANSF_TAG)
                elements['Beams'][self.sections['Storey'][ele]].append(eleid)

        print('[SUCCESS] Elements have been defined')
        return elements, base_cols

    def create_joints(self):
        """
        Creates Joint elements
        :return: None
        """
        joint_id = 0
        for st in range(self.g.nst):
            for bay in range(self.g.nbays + 1):
                eleid = int(f"{st + 1}{bay + 1}")
                nd1 = int(f"{st + 2}{bay + 1}10")
                nd2 = int(f"{st + 2}{bay + 1}20")
                nd3 = int(f"{st + 2}{bay + 1}30")
                nd4 = int(f"{st + 2}{bay + 1}40")
                ndc = int(str(nd1)[:-1])
                mat1 = 100000 + nd1
                mat2 = 100000 + nd2
                mat3 = 100000 + nd3
                mat4 = 100000 + nd4
                matc = 200000 + joint_id
                if st == self.g.nst - 1:
                    mat3 = 0
                if bay == 0:
                    mat4 = 0
                if bay == self.g.nbays:
                    mat2 = 0
                op.element('Joint2D', eleid, nd1, nd2, nd3, nd4, ndc, mat1, mat2, mat3, mat4, matc, 1)
                joint_id += 1

        print('[SUCCESS] Element connectivity and joint elements have been defined')

    def define_pdelta_columns(self, option='Truss'):
        """
        Defines pdelta columns
        :param option: str                              Option for linking the gravity columns (Truss or EqualDOF)
        :return: None
        """
        option = option.lower()
        # Elastic modulus of concrete
        young_modulus = float(self.materials['Ec']) * 1000.0
        # Check whether Pdelta forces were provided (if not, skips step)
        if 'pdelta' in list(self.loads['Pattern']):
            # Material definition
            pdelta_mat_tag = 300000 if self.hingeModel == 'haselton' else self.g.nbays + 2
            if self.system == 'Perimeter':
                op.uniaxialMaterial('Elastic', pdelta_mat_tag, young_modulus)

            # X coordinate of the columns
            x_coord = self.g.widths[(-1)] + 3.0

            # Geometric transformation for the columns
            pdelta_transf_tag = 3
            op.geomTransf('Linear', pdelta_transf_tag)

            # Node creations and linking to the lateral load resisting structure
            for st in range(self.g.nst + 1):
                if st == 0:
                    if self.hingeModel == 'haselton':
                        node = int(f"1{self.g.nbays + 2}000")
                    else:
                        node = int(f"{pdelta_mat_tag}{st}")
                    # Create and fix the node
                    op.node(node, x_coord, self.g.heights[st])
                    op.fix(node, 1, 1, 0)
                else:
                    if self.hingeModel == 'haselton':
                        nodeFrame = int(f"{st + 1}{self.g.nbays + 1}20")
                        node = int(f"{st}{self.g.nbays + 2}")
                        ele = int(f"2{st}{self.g.nbays + 1}")
                    else:
                        nodeFrame = int(f"{self.g.nbays + 1}{st}")
                        node = int(f"{pdelta_mat_tag}{st}")
                        ele = int(f"1{self.g.nbays + 1}{st}")
                    # Create the node
                    op.node(node, x_coord, self.g.heights[st])

                    if option == 'truss':
                        op.element('Truss', ele, nodeFrame, node, 5.0, pdelta_mat_tag)

                    elif option == 'equaldof':
                        for bay in range(self.g.nbays + 1):
                            if self.hingeModel == 'haselton':
                                op.equalDOF(node, int(f"{st + 1}{bay + 1}1"), 1)
                            else:
                                op.equalDOF(node, int(f"{bay + 1}{st}"), 1)

                    else:
                        raise ValueError('[EXCEPTION] Wrong option for linking gravity columns (needs to be Truss '
                                         'or EqualDOF')

            # Creation of P-Delta column elements
            agcol = 0.5**2
            izgcol = 0.5**4/12/1.0e4
            for st in range(1, self.g.nst + 1):
                if self.hingeModel == 'haselton':
                    eleid = int(f"1{st}{self.g.nbays + 2}")
                    if st == 0:
                        node_i = int(f"1{self.g.nbays + 2}000")
                    else:
                        node_i = int(f"{st - 1}{self.g.nbays + 2}")
                    node_j = int(f"{st}{self.g.nbays + 2}")
                else:
                    eleid = int(f"2{pdelta_mat_tag}{st}")
                    node_i = int(f"{pdelta_mat_tag}{st - 1}")
                    node_j = int(f"{pdelta_mat_tag}{st}")
                op.element('elasticBeamColumn', eleid, node_i, node_j, agcol, young_modulus, izgcol, pdelta_transf_tag)

            # Definition of loads
            op.timeSeries('Linear', 11)
            op.pattern('Plain', 11, 11)
            pdelta_loads = self.loads[(self.loads['Pattern'] == 'pdelta')].reset_index(drop=True)
            for st in range(1, self.g.nst + 1):
                load = pdelta_loads[(pdelta_loads['Storey'] == st)]['Load'].iloc[0]
                if not load:
                    pass
                else:
                    if self.hingeModel == 'haselton':
                        op.load(int(f"{st}{self.g.nbays + 2}"), self.NEGLIGIBLE, -load, self.NEGLIGIBLE)
                    else:
                        op.load(int(f"{pdelta_mat_tag}{st}"), self.NEGLIGIBLE, -load, self.NEGLIGIBLE)

            print('[SUCCESS] P-Delta columns have been defined')

    def define_masses(self):
        """
        Defining masses
        :return: None
        """
        masses = self.loads[(self.loads['Pattern'] == 'mass')].reset_index(drop=True)
        for st in range(self.g.nst):
            for bay in range(self.g.nbays + 1):
                if bay == 0 or bay == self.g.nbays:
                    m = masses[(masses['Storey'] == st + 1)]['Load'].iloc[0] / (2 * self.g.nbays)
                else:
                    m = masses[(masses['Storey'] == st + 1)]['Load'].iloc[0] / self.g.nbays

                # Assign the masses
                if self.hingeModel == 'haselton':
                    op.mass(int(f"{st + 2}{bay + 1}10"), m, self.NEGLIGIBLE, self.NEGLIGIBLE)
                else:
                    op.mass(int(f"{bay + 1}{st + 1}"), m, self.NEGLIGIBLE, self.NEGLIGIBLE)

        else:
            print('[SUCCESS] Seismic masses have been defined')

    def set_recorders(self, analysis, **kwargs):
        """
        Defining recorders
        :param analysis: str
        :param kwargs:
        :return: dict                               Dictionary containing all results
        """
        results = {}
        r = Recorders(self.g.nbays, self.g.nst, self.elements, self.hingeModel)
        base_nodes = kwargs.get('base_nodes', None)
        num_modes = kwargs.get('num_modes', None)

        if analysis == 'ST' or analysis == 'static' or analysis == 'gravity':
            results = r.st_recorder(base_nodes)

        elif analysis == 'PO' or analysis == 'pushover':
            push_directory = self.outputsDir / 'PO'
            self.createFolder(push_directory)
            r.po_recorder(base_nodes, push_directory)

        elif analysis == 'MA' or analysis == 'modal':
            results = r.ma_recorder(num_modes)

        elif analysis == 'ELF' or analysis == 'ELFM':
            results = r.st_recorder(base_nodes)

        else:
            results = None

        print('[SUCCESS] Recorders have been generated')
        return results

    def define_loads(self, elements, apply_loads=True):
        """
        Defines gravity loads provided by the user
        :param elements: dict                           Dictionary containing IDs of all elements
        :param apply_loads: bool                        Whether to apply loads or not
        :return: None
        """
        if apply_loads:
            op.timeSeries('Linear', 1)
            op.pattern('Plain', 1, 1)
            distributed = self.loads[(self.loads['Pattern'] == 'distributed')].reset_index(drop=True)
            if self.hingeModel == 'haselton':
                for idx in range(1, self.g.nst + 1):
                    ele_ids = elements['Beams'][idx]
                    load = distributed[(distributed['Storey'] == idx)]['Load'].iloc[0]
                    if not ele_ids:
                        pass
                    else:
                        for ele in ele_ids:
                            op.eleLoad('-ele', ele, '-type', '-beamUniform', -load)

            else:
                for ele in elements['Beams']:
                    load = distributed[(distributed['Storey'] == int(ele[(-1)]))]['Load'].iloc[0]
                    op.eleLoad('-ele', int(ele), '-type', '-beamUniform', -load)

            print('[SUCCESS] Gravity loads have been defined')

    def perform_analysis(self, elfm_filename=None, **kwargs):
        """
        Performs analysis
        :param elfm_filename: str                           File name of ELF containing load values (for ELF)
        :param kwargs: mode_shape: list                     Mode shape values from MA (for PO)
        :param kwargs: damping: float                       Damping (for MA)
        :param kwargs: spo_pattern: int                     SPO lateral load pattern ID
        :return: dict                                       Results containing recorded data
        """
        mode_shape = kwargs.get('mode_shape', None)
        damping = kwargs.get('damping', None)
        spo_pattern = kwargs.get('spo_pattern', 2)
        control_nodes = []
        for i in range(self.g.nst):
            if self.hingeModel == 'haselton':
                control_nodes.append(int(f"{i + 2}{self.g.nbays + 1}20"))
            else:
                control_nodes.append(int(f"{self.g.nbays + 1}{i + 1}"))

        if 'ELF' in self.analysis_type or 'ELFM' in self.analysis_type:
            print('[STEP] ELF started')
            try:
                elfm_forces = pd.read_csv(elfm_filename)
                elfm_forces = elfm_forces[(elfm_forces['Pattern'] == 'elf')]
            except:
                raise ValueError('[EXCEPTION] ELF forces not provided')

            op.timeSeries('Linear', 3)
            op.pattern('Plain', 300, 3)
            for st in range(self.g.nst):
                load = elfm_forces[(elfm_forces['Storey'] == st + 1)]['Load'].iloc[0]
                if self.hingeModel == 'haselton':
                    op.load(int(f"{st + 1}{self.g.nbays + 1}"), load, self.NEGLIGIBLE, self.NEGLIGIBLE)
                else:
                    op.load(int(f"{self.g.nbays + 1}{st + 1}"), load, self.NEGLIGIBLE, self.NEGLIGIBLE)
            else:
                s = Static()
                s.static_analysis()
                self.results['ELF'] = self.set_recorders('ELF', base_nodes=self.base_nodes)
                print('[SUCCESS] ELF done')

        if 'ST' in self.analysis_type or 'static' in self.analysis_type or 'gravity' in self.analysis_type:
            print('[STEP] Gravity static analysis started')
            s = Static()
            s.static_analysis()
            self.results['Gravity'] = self.set_recorders('ST', base_nodes=self.base_nodes)

            filepath = self.outputsDir / 'ST'

            with open(f"{filepath}.json", 'w') as (f):
                json.dump(self.results['Gravity'], f)

            print('[SUCCESS] Static gravity analysis done')

        if 'MA' in self.analysis_type or 'modal' in self.analysis_type:
            print('[STEP] Modal analysis started')
            m = Modal(self.NUM_MODES, self.DAMP_MODES, damping)
            self.results['Modal'] = self.set_recorders('MA', num_modes=self.NUM_MODES)
            self.results['Modal']['Periods'] = m.period
            self.results['Modal']['Damping'] = m.xi_modes
            self.results['Modal']['CircFreq'] = m.omega
            filepath = self.outputsDir / 'MA'
            with open(f"{filepath}.json", 'w') as (f):
                json.dump(self.results['Modal'], f)
            print('[SUCCESS] Modal analysis done')
            mode_shape = self.results['Modal']['Mode1']

        if 'PO' in self.analysis_type or 'pushover' in self.analysis_type:
            print('[STEP] Static pushover analysis started')
            if self.hingeModel == 'haselton':
                id_ctrl_node = int(f"{self.g.nst + 1}{self.g.nbays + 1}20")
            else:
                id_ctrl_node = int(f"{self.g.nbays + 1}{self.g.nst}")
            id_ctrl_dof = 1
            spo = SPO(id_ctrl_node, id_ctrl_dof, self.base_cols)
            spo.load_pattern(control_nodes, load_pattern=spo_pattern, heights=self.g.heights, mode_shape=mode_shape)
            spo.set_analysis()
            outputs = spo.seek_solution()
            filepath = self.outputsDir / 'SPO'
            with open(f"{filepath}.pickle", 'wb') as (f):
                pickle.dump(outputs, f)
            print('[SUCCESS] Static pushover analysis done')

    def lumped_hinge_element(self):
        """
        Creates lumped hinge elements for the hysteretic model
        :return: dict                       Dictionary containing element IDs to be used for recording internal forces
        """
        s = Sections(self.sections, self.materials)
        elements = {'Columns': [], 'Beams': []}
        base_cols = []
        for index, ele in self.sections.iterrows():
            if ele['Element'] == 'Beam':
                eleid = f"1{ele['Bay']}{ele['Storey']}"
                transfTag = self.BEAM_TRANSF_TAG
                elements['Beams'].append(eleid)
            else:
                eleid = f"2{ele['Bay']}{ele['Storey']}"
                transfTag = self.COL_TRANSF_TAG
                elements['Columns'].append(eleid)
                if ele['Storey'] == 1:
                    base_cols.append(eleid)
            s.hysteretic_hinges(ele, transfTag)

        return elements, base_cols

    def model(self):
        """
        Creates the full model
        :return: None
        """
        self.create_model()
        self.base_nodes, hinge_nodes = self.create_nodes()
        self.define_transformations()
        if self.hingeModel == 'haselton':
            self.joint_materials()
            self.rot_springs(self.base_nodes)
            self.bilin_springs(hinge_nodes)
            self.elements, self.base_cols = self.create_elements()
            self.create_joints()
        elif self.hingeModel == 'hysteretic':
            self.elements, self.base_cols = self.lumped_hinge_element()
        else:
            raise ValueError('[EXCEPTION] Wrong lumped hinge model (should be Hysteretic or Haselton)')

        self.define_pdelta_columns(option='EqualDOF')
        self.define_masses()
