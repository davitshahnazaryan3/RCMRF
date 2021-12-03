"""
Model creator of an RC MRF. Lumped hinge models following the recommendations of Haselton 2007 are used.
"""
import openseespy.opensees as op
import pandas as pd

from client.geometry import Geometry
from client.sections import Sections
from client.recorders import Recorders
from analysis.static import Static
from analysis.modal import Modal
from analysis.spo import SPO
from utils.utils import *


class ModelToTCL:
    def __init__(self, analysis_type, sections_file, loads_file, materials, outputsDir, system='space',
                 hingeModel='hysteretic', flag3d=False, direction=0, tcl_filename=None):
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
        :param flag3d: bool                         True for 3D modelling, False for 2D modelling
        :param direction: int                       Direction of application, 0: x; 1: y
        :param tcl_filename: str                    TCL filename to export model to
        """
        self.base_nodes = None
        self.base_cols = None
        self.elements = None
        self.analysis_type = analysis_type
        self.materials = pd.read_csv(materials)
        self.system = system.lower()
        self.hingeModel = hingeModel.lower()
        self.flag3d = flag3d
        self.direction = direction
        self.COL_TRANSF_TAG = 1
        self.BEAM_X_TRANSF_TAG = 2
        self.BEAM_Y_TRANSF_TAG = 3
        self.NEGLIGIBLE = 1.e-09
        self.UBIG = 1.e10
        self.outputsDir = outputsDir

        # Nodes necessary for SPO analysis
        self.spo_nodes = []
        if self.flag3d:
            f = {}
            cols = None
            for i in sections_file:
                try:
                    f[i] = pd.read_csv(sections_file[i])
                except:
                    f[i] = sections_file[i]
                cols = [j for j in f[i].columns if j not in ('Element', 'Position', 'Storey', 'Bay', "Direction")]

                for col in cols:
                    f[i][col] = f[i][col].astype(float)

            # Rename keys of f
            try:
                f["x"] = f.pop("x_seismic")
                f["y"] = f.pop("y_seismic")
            except:
                pass

            self.sections = f

        else:
            try:
                self.sections = pd.read_csv(sections_file)
            except:
                self.sections = sections_file

            cols = [i for i in self.sections.columns if i not in ('Element', 'Position', 'Storey', 'Bay', "Direction")]

            for col in cols:
                self.sections[col] = self.sections[col].astype(float)

        self.loads = pd.read_csv(loads_file)
        check_integrity(self.system, self.flag3d, self.hingeModel)
        self.g = Geometry(self.sections, self.hingeModel, flag3d=self.flag3d)
        self.NUM_MODES = 3
        self.DAMP_MODES = [1, 2]
        self.results = {}

        # File to export model to
        self.file = None
        self.tcl_filename = tcl_filename

    def create_model(self):
        """
        Initiates model creation
        :return: None
        """
        createFolder(self.outputsDir / "Models")

        if self.flag3d:
            op.model('Basic', '-ndm', 3, '-ndf', 6)

            # Write to tcl file
            if "PO" in self.analysis_type:
                d = "x" if self.direction == 0 else "y"
                filename = self.outputsDir / f"Models/{self.tcl_filename}_pushover_{d}.tcl"
            elif 'ST' in self.analysis_type or 'static' in self.analysis_type or 'gravity' in self.analysis_type:
                filename = self.outputsDir / f"Models/{self.tcl_filename}_static.tcl"
            elif 'MA' in self.analysis_type or 'modal' in self.analysis_type:
                filename = self.outputsDir / f"Models/{self.tcl_filename}_modal.tcl"
            else:
                filename = self.outputsDir / f"Models/{self.tcl_filename}.tcl"

            lines = ["# Create Model Global", "wipe;", "model BasicBuilder -ndm 3 -ndf 6;"]
            self.file = open(filename, "w+")
            self.file.write("\n".join(lines))

        else:
            op.model('Basic', '-ndm', 2, '-ndf', 3)
        print('[INITIATE] Model generation started')

    def create_nodes(self, fixity='fixed'):
        """
        Creates nodes
        :param fixity: str                          Boundary condition of base nodes
        :return: list                               Base node IDs
        """
        # Fixity for ground floor columns of seismic frames
        if fixity == 'fixed':
            fix = 1
        elif fixity == 'pinned':
            fix = 0
        else:
            raise Exception('[EXCEPTION] Wrong boundary condition for base nodes')

        # Number of bays in x and y directions, spans
        if self.flag3d:
            spans_x = self.g.widths[0]
            spans_y = self.g.widths[1]

        else:
            spans_x = spans_y = self.g.widths

        # Get dataframe containing all nodes of the building / frame
        df = self.g.define_nodes()

        self.file.write("\n\n# Define nodes")
        # Initialization of base nodes and hinges
        base_nodes = []
        hinge_nodes = []
        for n in df.index:
            nodetag = int(df['Node id'][n])
            xloc = df['x'][n]
            zloc = df['z'][n]

            # Defining nodes
            if self.flag3d:
                yloc = df['y'][n]
                op.node(nodetag, xloc, yloc, zloc)
                self.file.write(f"\nnode {nodetag} {xloc} {yloc} {zloc};")
            else:
                yloc = 0.
                op.node(int(df['Node id'][n]), xloc, zloc)

            # Restraining ground floor nodes
            if zloc == 0 and nodetag < 10000:
                if self.flag3d:
                    if self.system == "perimeter":
                        # Fix or pin the base nodes
                        if (xloc == 0. or xloc == max(spans_x)) and (yloc == 0. or yloc == max(spans_y)):
                            # Fix the external columns of seismic frames
                            op.fix(nodetag, 1, 1, 1, fix, fix, fix)
                        elif 0. < xloc < max(spans_x) and (yloc == 0. or yloc == max(spans_y)):
                            # Fix internal columns of X seismic frames against rotation in x direction
                            op.fix(nodetag, 1, 1, 1, 0, fix, 0)
                        elif 0. < yloc < max(spans_y) and (xloc == 0. or xloc == max(spans_x)):
                            # Fix internal columns of Y seismic frames against rotation in y direction
                            op.fix(nodetag, 1, 1, 1, fix, 0, 0)
                        else:
                            # Pin the columns of gravity columns
                            op.fix(nodetag, 1, 1, 1, 0, 0, 0)
                    else:
                        # Fix or pin all base nodes
                        op.fix(nodetag, 1, 1, 1, fix, fix, fix)
                        self.file.write(f"\nfix {nodetag} 1 1 1 1 1 1;")

                else:
                    # Fix the base columns
                    op.fix(nodetag, 1, 1, fix)

                base_nodes.append(nodetag)

            elif 0.0 <= xloc <= max(spans_x):
                # Not necessary for hysteretic hinges
                if nodetag < 10000 and zloc <= max(self.g.heights):
                    hinge_nodes.append(nodetag)

        print('[SUCCESS] Geometric properties have been defined')
        return base_nodes, hinge_nodes

    def define_transformations(self, col_transf_type='PDelta', beam_transf_tag='PDelta'):
        """
        Defines geometric transformations for beams and columns (PDelta, Linear, Corotational)
        :param col_transf_type: str                 Column transformation type
        :param beam_transf_tag: str                 Beam transformation type
        :return: None
        """
        if any((tag not in ('PDelta', 'Linear', 'Corotational') for tag in (col_transf_type, beam_transf_tag))):
            raise Exception('[EXCEPTION] Wrong transformation type provided')

        else:
            if self.flag3d:
                # TODO add logic for precise estimation of offsets
                op.geomTransf(col_transf_type, self.COL_TRANSF_TAG, 0., -1., 0.,
                              "-jntOffset", 0.0, 0.0, 0.3, 0.0, 0.0, -0.3)
                op.geomTransf(beam_transf_tag, self.BEAM_X_TRANSF_TAG, 0., -1., 0.,
                              "-jntOffset", 0.2, 0.0, 0.0, -0.2, 0.0, 0.0)
                op.geomTransf(beam_transf_tag, self.BEAM_Y_TRANSF_TAG, 1., 0., 0.,
                              "-jntOffset", 0.0, 0.2, 0.0, 0.0, -0.2, 0.0)

                self.file.write("\n\n# Define geometric transformations")
                self.file.write(f"\ngeomTransf {col_transf_type} {self.COL_TRANSF_TAG} 0. -1. 0."
                                f" -jntOffset 0.0 0.0 0.3 0.0 0.0 -0.3;")
                self.file.write(f"\ngeomTransf {beam_transf_tag} {self.BEAM_X_TRANSF_TAG} 0. -1. 0."
                                f" -jntOffset 0.2 0.0 0.0 -0.2 0.0 0.0;")
                self.file.write(f"\ngeomTransf {beam_transf_tag} {self.BEAM_Y_TRANSF_TAG} 1. 0. 0."
                                f" -jntOffset 0.0 0.2 0.0 0.0 -0.2 0.0;")

            else:
                op.geomTransf(col_transf_type, self.COL_TRANSF_TAG)
                op.geomTransf(beam_transf_tag, self.BEAM_X_TRANSF_TAG)

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
                op.element('elasticBeamColumn', eleid, node_i, node_j, area, young_modulus, inertia,
                           self.COL_TRANSF_TAG)
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
                op.element('elasticBeamColumn', eleid, node_i, node_j, area, young_modulus, inertia,
                           self.COL_TRANSF_TAG)
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
                op.element('elasticBeamColumn', eleid, node_i, node_j, area, young_modulus, inertia,
                           self.BEAM_X_TRANSF_TAG)
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
            pdelta_mat_tag = 300000 if self.hingeModel == 'haselton' else int(self.g.nbays + 2)
            if self.system == 'perimeter':
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
        self.file.write("\n\n# Define masses")

        if self.flag3d:
            nbays_x = max(self.sections["x"]["Bay"] - 1)
            nbays_y = max(self.sections["y"]["Bay"] - 1)
            nst = max(self.sections["x"]["Storey"])
            spans_x = np.diff(self.g.widths[0])
            spans_y = np.diff(self.g.widths[1])
            for st in range(1, nst + 1):
                for xbay in range(1, nbays_x + 2):
                    for ybay in range(1, nbays_y + 2):
                        nodetag = int(f"{xbay}{ybay}{st}")

                        # Corner nodes
                        if xbay == 1:
                            if ybay == 1:
                                # Corner node
                                area = spans_x[xbay - 1] * spans_y[ybay - 1] / 4
                            elif ybay == nbays_y + 1:
                                # Corner node
                                area = spans_x[xbay - 1] * spans_y[ybay - 2] / 4
                            else:
                                # Side node
                                area = spans_x[xbay - 1] * (spans_y[ybay - 2] + spans_y[ybay - 1]) / 4

                        elif xbay == nbays_x + 1:
                            if ybay == 1:
                                # Corner node
                                area = spans_x[xbay - 2] * spans_y[ybay - 1] / 4
                            elif ybay == nbays_y + 1:
                                # Corner node
                                area = spans_x[xbay - 2] * spans_y[ybay - 2] / 4
                            else:
                                # Side node
                                area = spans_x[xbay - 2] * (spans_y[ybay - 2] + spans_y[ybay - 1]) / 4

                        else:
                            if ybay == 1:
                                # Side node
                                area = (spans_x[xbay - 2] + spans_x[xbay - 1]) * spans_y[ybay - 1] / 4
                            elif ybay == nbays_y + 1:
                                # Side node
                                area = (spans_x[xbay - 2] + spans_x[xbay - 1]) * spans_y[ybay - 2] / 4
                            else:
                                # Internal node
                                area = (spans_x[xbay - 2] + spans_x[xbay - 1]) * (
                                            spans_y[ybay - 2] + spans_y[ybay - 1]) / 4

                        # Mass based on tributary area
                        q = self.loads[(self.loads["Pattern"] == "seismic") &
                                       (self.loads["Storey"] == st)]["Load"].iloc[0]
                        mass = area * q / 9.81
                        op.mass(nodetag, mass, mass, mass, self.NEGLIGIBLE, self.NEGLIGIBLE, self.NEGLIGIBLE)
                        self.file.write(f"\nmass {nodetag} {mass} {mass} {mass} {self.NEGLIGIBLE} {self.NEGLIGIBLE} "
                                        f"{self.NEGLIGIBLE};")

        else:
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
        r = Recorders(self.g, self.elements, self.hingeModel, self.flag3d)
        base_nodes = kwargs.get('base_nodes', None)
        num_modes = kwargs.get('num_modes', None)

        if analysis == 'ST' or analysis == 'static' or analysis == 'gravity':
            results = r.st_recorder(base_nodes)

        elif analysis == 'MA' or analysis == 'modal':
            lam = kwargs.get('lam', None)
            results = r.ma_recorder(num_modes, lam, self.outputsDir)

        elif analysis == 'ELF' or analysis == 'ELFM':
            results = r.st_recorder(base_nodes)

        else:
            results = None

        print('[SUCCESS] Recorders have been generated')
        return results

    def define_loads(self, elements, apply_loads=True, apply_point=False):
        """
        Defines gravity loads provided by the user
        :param elements: dict                           Dictionary containing IDs of all elements
        :param apply_loads: bool                        Whether to apply loads or not
        :param apply_point: bool                        Whether to apply loads as point loads (if both True, advantage
                                                        will be given to point loads)
        :return: None
        """
        self.file.write("\n\n# Apply gravity loads")

        # For now, point loads are not created for Haselton model, so force distributed loads
        if self.hingeModel == "haselton":
            apply_loads = True
            apply_point = False

        if apply_loads:
            op.timeSeries('Linear', 1)
            op.pattern('Plain', 1, 1)
            self.file.write("\npattern Plain 1 Linear {")

            if self.hingeModel == 'haselton':
                distributed = self.loads[(self.loads['Pattern'] == 'distributed')].reset_index(drop=True)
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
                    if self.flag3d:
                        loads = self.loads[(self.loads['Pattern'] == 'q')].reset_index(drop=True)
                        spans_x = np.diff(self.g.widths[0])
                        spans_y = np.diff(self.g.widths[1])

                        for beam in elements["Beams"][ele]:
                            st = int(str(beam)[-1])
                            xbay = int(str(beam)[1])
                            ybay = int(str(beam)[2])
                            q = loads[loads["Storey"] == st]["Load"].iloc[0]

                            if ele == "x" or ele == "gravity_x":
                                # Beams along X direction
                                # Load over a beam
                                control_length = spans_y[ybay - 1] if ybay < len(spans_y) + 1 else spans_y[ybay - 2]
                                if spans_x[xbay - 1] <= control_length:
                                    # Triangular rule
                                    load = q * spans_x[xbay - 1] ** 2 / 4 / spans_x[xbay - 1]
                                else:
                                    # Trapezoidal rule
                                    load = 1 / 4 * q * spans_y[ybay - 1] * (
                                                2 * spans_x[xbay - 1] - spans_y[ybay - 1]) / spans_x[xbay - 1]
                                load = round(load, 2)

                                # End nodes
                                nodei = beam - 3000
                                nodej = beam - 3000 + 100

                                if apply_point:
                                    op.load(nodei, self.NEGLIGIBLE, self.NEGLIGIBLE, -load * spans_x[xbay - 1] / 2,
                                            self.NEGLIGIBLE, self.NEGLIGIBLE, self.NEGLIGIBLE)
                                    op.load(nodej, self.NEGLIGIBLE, self.NEGLIGIBLE, -load * spans_x[xbay - 1] / 2,
                                            self.NEGLIGIBLE, self.NEGLIGIBLE, self.NEGLIGIBLE)
                                    self.file.write(f"\n\tload {nodei} {self.NEGLIGIBLE} {self.NEGLIGIBLE}"
                                                    f" {-load * spans_x[xbay - 1] / 2} {self.NEGLIGIBLE}"
                                                    f" {self.NEGLIGIBLE} {self.NEGLIGIBLE};")
                                    self.file.write(f"\n\tload {nodej} {self.NEGLIGIBLE} {self.NEGLIGIBLE}"
                                                    f" {-load * spans_x[xbay - 1] / 2} {self.NEGLIGIBLE}"
                                                    f" {self.NEGLIGIBLE} {self.NEGLIGIBLE};")

                                else:
                                    op.eleLoad('-ele', beam, '-type', '-beamUniform', -load, self.NEGLIGIBLE)
                                    self.file.write(f"\n\teleLoad -ele {beam} -type -beamUniform -{load} "
                                                    f"{self.NEGLIGIBLE};")

                                # Additional load for interior beams
                                if 1 < ybay < len(spans_y) + 1:
                                    if spans_x[xbay - 1] <= spans_y[ybay - 2]:
                                        # Triangular rule
                                        load = q * spans_x[xbay - 1] ** 2 / 4 / spans_x[xbay - 1]
                                    else:
                                        # Trapezoidal rule
                                        load = 1 / 4 * q * spans_y[ybay - 2] * (
                                                    2 * spans_x[xbay - 1] - spans_y[ybay - 2]) / spans_x[xbay - 1]
                                    load = round(load, 2)

                                    # Applying the load
                                    if apply_point:
                                        op.load(nodei, self.NEGLIGIBLE, self.NEGLIGIBLE,
                                                -load * spans_x[xbay - 1] / 2,
                                                self.NEGLIGIBLE, self.NEGLIGIBLE, self.NEGLIGIBLE)
                                        op.load(nodej, self.NEGLIGIBLE, self.NEGLIGIBLE,
                                                -load * spans_x[xbay - 1] / 2,
                                                self.NEGLIGIBLE, self.NEGLIGIBLE, self.NEGLIGIBLE)
                                        self.file.write(f"\n\tload {nodei} {self.NEGLIGIBLE} {self.NEGLIGIBLE}"
                                                        f" {-load * spans_x[xbay - 1] / 2} {self.NEGLIGIBLE}"
                                                        f" {self.NEGLIGIBLE} {self.NEGLIGIBLE};")
                                        self.file.write(f"\n\tload {nodej} {self.NEGLIGIBLE} {self.NEGLIGIBLE}"
                                                        f" {-load * spans_x[xbay - 1] / 2} {self.NEGLIGIBLE}"
                                                        f" {self.NEGLIGIBLE} {self.NEGLIGIBLE};")

                                    else:
                                        op.eleLoad('-ele', beam, '-type', '-beamUniform', -load, self.NEGLIGIBLE)
                                        self.file.write(f"\n\teleLoad -ele {beam} -type -beamUniform -{load} "
                                                        f"{self.NEGLIGIBLE};")

                            else:
                                # Beams along Y direction
                                # Load over a beam
                                control_length = spans_x[xbay - 1] if xbay < len(spans_x) + 1 else spans_x[xbay - 2]
                                if spans_y[ybay - 1] <= control_length:
                                    # Triangular rule
                                    load = q * spans_y[ybay - 1] ** 2 / 4 / spans_y[ybay - 1]
                                else:
                                    # Trapezoidal rule
                                    load = 1 / 4 * q * spans_x[xbay - 1] * \
                                           (2 * spans_y[ybay - 1] - spans_x[xbay - 1]) / spans_y[ybay - 1]
                                load = round(load, 2)

                                # End nodes
                                nodei = beam - 2000
                                nodej = beam - 2000 + 10

                                # Applying the load
                                if apply_point:
                                    op.load(nodei, self.NEGLIGIBLE, self.NEGLIGIBLE, -load * spans_y[ybay - 1] / 2,
                                            self.NEGLIGIBLE, self.NEGLIGIBLE, self.NEGLIGIBLE)
                                    op.load(nodej, self.NEGLIGIBLE, self.NEGLIGIBLE, -load * spans_y[ybay - 1] / 2,
                                            self.NEGLIGIBLE, self.NEGLIGIBLE, self.NEGLIGIBLE)
                                    self.file.write(f"\n\tload {nodei} {self.NEGLIGIBLE} {self.NEGLIGIBLE}"
                                                    f" {-load * spans_y[ybay - 1] / 2} {self.NEGLIGIBLE}"
                                                    f" {self.NEGLIGIBLE} {self.NEGLIGIBLE};")
                                    self.file.write(f"\n\tload {nodej} {self.NEGLIGIBLE} {self.NEGLIGIBLE}"
                                                    f" {-load * spans_y[ybay - 1] / 2} {self.NEGLIGIBLE}"
                                                    f" {self.NEGLIGIBLE} {self.NEGLIGIBLE};")

                                else:
                                    op.eleLoad('-ele', beam, '-type', '-beamUniform', -load, self.NEGLIGIBLE)
                                    self.file.write(f"\n\teleLoad -ele {beam} -type -beamUniform -{load} "
                                                    f"{self.NEGLIGIBLE};")

                                # Additional load for interior beams
                                if 1 < xbay < len(spans_x) + 1:
                                    if spans_y[ybay - 1] <= spans_x[xbay - 2]:
                                        # Triangular rule
                                        load = q * spans_y[ybay - 1] ** 2 / 4 / spans_y[ybay - 1]
                                    else:
                                        # Trapezoidal rule
                                        load = 1 / 4 * q * spans_x[xbay - 2] * \
                                               (2 * spans_y[ybay - 1] - spans_x[xbay - 2]) / spans_y[ybay - 1]
                                    load = round(load, 2)

                                    # Applying the load
                                    if apply_point:
                                        op.load(nodei, self.NEGLIGIBLE, self.NEGLIGIBLE,
                                                -load * spans_y[ybay - 1] / 2,
                                                self.NEGLIGIBLE, self.NEGLIGIBLE, self.NEGLIGIBLE)
                                        op.load(nodej, self.NEGLIGIBLE, self.NEGLIGIBLE,
                                                -load * spans_y[ybay - 1] / 2,
                                                self.NEGLIGIBLE, self.NEGLIGIBLE, self.NEGLIGIBLE)
                                        self.file.write(f"\n\tload {nodei} {self.NEGLIGIBLE} {self.NEGLIGIBLE}"
                                                        f" {-load * spans_y[ybay - 1] / 2} {self.NEGLIGIBLE}"
                                                        f" {self.NEGLIGIBLE} {self.NEGLIGIBLE};")
                                        self.file.write(f"\n\tload {nodej} {self.NEGLIGIBLE} {self.NEGLIGIBLE}"
                                                        f" {-load * spans_y[ybay - 1] / 2} {self.NEGLIGIBLE}"
                                                        f" {self.NEGLIGIBLE} {self.NEGLIGIBLE};")
                                    else:
                                        op.eleLoad('-ele', beam, '-type', '-beamUniform', -load, self.NEGLIGIBLE)
                                        self.file.write(f"\n\teleLoad -ele {beam} -type -beamUniform -{load} "
                                                        f"{self.NEGLIGIBLE};")

                    else:
                        distributed = self.loads[(self.loads['Pattern'] == 'distributed')].reset_index(drop=True)
                        load = distributed[(distributed['Storey'] == int(ele[-1]))]['Load'].iloc[0]
                        if apply_point:
                            # Storey and bay levels
                            st = int(ele[-1])
                            bay = int(ele[1])
                            # Bay width
                            w = self.g.beams[(self.g.beams['Bay'] == bay)].iloc[0]['length']
                            # Point load value
                            p = load / 2 * w
                            # Apply the loads
                            op.load(int(f"{bay}{st}"), self.NEGLIGIBLE, -p, self.NEGLIGIBLE)
                            op.load(int(f"{bay+1}{st}"), self.NEGLIGIBLE, -p, self.NEGLIGIBLE)

                        else:
                            # Distributed load (kN/m)
                            op.eleLoad('-ele', int(ele), '-type', '-beamUniform', -load)

            if apply_point:
                print('[SUCCESS] Gravity loads as point loads have been defined')
            else:
                print('[SUCCESS] Gravity loads aas distributed loads have been defined')

            self.file.write("\n};")

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
                    if self.flag3d:
                        # TODO, to be corrected for 3D
                        op.load(int(f"{self.g.nbays + 1}{st + 1}"), load, self.NEGLIGIBLE, self.NEGLIGIBLE,
                                self.NEGLIGIBLE, self.NEGLIGIBLE, self.NEGLIGIBLE)
                    else:
                        op.load(int(f"{self.g.nbays + 1}{st + 1}"), load, self.NEGLIGIBLE, self.NEGLIGIBLE)

            s = Static()
            s.static_analysis(self.flag3d)
            self.results['ELF'] = self.set_recorders('ELF', base_nodes=self.base_nodes)
            print('[SUCCESS] ELF done')

        if 'ST' in self.analysis_type or 'static' in self.analysis_type or 'gravity' in self.analysis_type:
            print('[STEP] Gravity static analysis started')
            s = Static()
            s.static_analysis(self.outputsDir, self.flag3d)
            self.results['Gravity'] = self.set_recorders('ST', base_nodes=self.base_nodes)

            filepath = self.outputsDir / 'ST'

            with open(f"{filepath}.json", 'w') as (f):
                json.dump(self.results['Gravity'], f)

            print('[SUCCESS] Static gravity analysis done')

        if 'MA' in self.analysis_type or 'modal' in self.analysis_type:
            self.file.write("\n\n# Call Modal analysis")
            self.file.write("\nsource modal_recorders.tcl")
            self.file.write("\nsource modal_analysis.tcl")
            self.file.write("\nwipe;")

            print('[STEP] Modal analysis started')
            m = Modal(self.NUM_MODES, self.DAMP_MODES, damping, self.outputsDir)
            self.results['Modal'], positions = self.set_recorders('MA', num_modes=self.NUM_MODES, lam=m.lam)

            # Modify positions of modal parameters
            if positions is not None:
                omega = [m.omega[i] for i in positions]
                xi_modes = [m.xi_modes[i] for i in positions]
                period = [m.period[i] for i in positions]
            else:
                omega = m.omega
                xi_modes = m.xi_modes
                period = m.period

            self.results['Modal']['Periods'] = period
            self.results['Modal']['Damping'] = xi_modes
            self.results['Modal']['CircFreq'] = omega

            filepath = self.outputsDir / 'MA'
            with open(f"{filepath}.json", 'w') as (f):
                json.dump(self.results['Modal'], f)

            print('[SUCCESS] Modal analysis done')
            mode_shape = self.results['Modal']['Mode1']

        if 'PO' in self.analysis_type or 'pushover' in self.analysis_type:
            self.file.write("\n\n# Static analysis")
            self.file.write("\nsource static.tcl")
            self.file.write("\n\n# Call Pushover analysis")
            d = "x" if self.direction == 0 else "y"
            self.file.write(f"\nsource spo_recorders_{d}_{self.tcl_filename[6:]}.tcl")
            self.file.write(f"\nsource spo_analysis_{d}_{self.tcl_filename[6:]}.tcl")
            self.file.write("\nwipe;")

            control_nodes = []
            for i in range(self.g.nst):
                if not self.flag3d:
                    if self.hingeModel == 'haselton':
                        control_nodes.append(int(f"{i + 2}{self.g.nbays + 1}20"))
                    else:
                        control_nodes.append(int(f"{self.g.nbays + 1}{i + 1}"))
                else:
                    control_nodes.append(int(f"{self.g.nbays[0] + 1}{self.g.nbays[1] + 1}{i + 1}"))

            print('[STEP] Static pushover analysis started')
            # Top node tag for recording the top displacement and for the integrator
            if not self.flag3d:
                if self.hingeModel == 'haselton':
                    id_ctrl_node = int(f"{self.g.nst + 1}{self.g.nbays + 1}20")
                else:
                    id_ctrl_node = int(f"{self.g.nbays + 1}{self.g.nst}")
                # Number of bays
                nbays_x = self.g.nbays
                nbays_y = None
            else:
                id_ctrl_node = int(f"{self.g.nbays[0] + 1}{self.g.nbays[1] + 1}{self.g.nst}")
                nbays_x = self.g.nbays[0]
                nbays_y = self.g.nbays[1]

            # DOF associated with the direction of interest
            id_ctrl_dof = self.direction + 1

            # Reference displacement to which cycles are run
            dref = 0.1 * max(self.g.heights)

            # Call the SPO object
            spo = SPO(id_ctrl_node, id_ctrl_dof, self.base_cols, self.base_nodes, dref, flag3d=self.flag3d,
                      direction=self.direction, filename=self.outputsDir / "Models", site=self.tcl_filename[6:])

            spo.load_pattern(control_nodes, load_pattern=spo_pattern, heights=self.g.heights, mode_shape=mode_shape,
                             nbays_x=nbays_x, nbays_y=nbays_y)
            spo.set_analysis(heights=self.g.heights)
            outputs = spo.seek_solution()
            filepath = self.outputsDir / f'SPO_{self.tcl_filename[6:]}_{self.direction+1}'
            with open(f"{filepath}.pickle", 'wb') as (f):
                pickle.dump(outputs, f)

            print('[SUCCESS] Static pushover analysis done')

    def lumped_hinge_element(self):
        """
        Creates lumped hinge elements for the hysteretic model
        :return: dict                       Dictionary containing element IDs to be used for recording internal forces
        """
        s = Sections(self.sections, self.materials)
        self.file.write("\n\n# Define nonlinear elements")

        base_cols = []
        if self.flag3d:
            # Initialize elements
            elements = {'Columns': {"x": [], "y": [], "gravity": []},
                        'Beams': {"x": [], "y": [], "gravity_x": [], "gravity_y": []}}
            # Number of bays and storeys
            nbays_x = max(self.sections["x"]["Bay"] - 1)
            nbays_y = max(self.sections["y"]["Bay"] - 1)
            nst = max(self.sections["x"]["Storey"])
            hinge_x = self.sections["x"]
            hinge_y = self.sections["y"]
            hinge_gr = self.sections["gravity"]

            # Add column elements
            for xbay in range(1, int(nbays_x + 2)):
                for ybay in range(1, int(nbays_y + 2)):
                    for st in range(1, int(nst + 1)):
                        previous_st = st - 1
                        # Element tag
                        et = int(f"1{xbay}{ybay}{st}")
                        # End nodes of column
                        inode = int(f"{xbay}{ybay}{previous_st}")
                        jnode = int(f"{xbay}{ybay}{st}")

                        # Transformation tag
                        transfTag = self.COL_TRANSF_TAG

                        # Columns of seismic frame along x direction
                        if ybay == 1 or ybay == nbays_y + 1:
                            # Columns of seismic frames along x direction
                            eleHinge = hinge_x[(hinge_x["Element"] == "Column") & (hinge_x["Bay"] == xbay) & (
                                    hinge_x["Storey"] == st)].reset_index(drop=True).iloc[0]
                            elements["Columns"]["x"].append(et)

                        elif (xbay == 1 or xbay == nbays_x + 1) and (1 < ybay < nbays_y + 1):
                            # Columns of seismic frames along y direction
                            eleHinge = hinge_y[(hinge_y["Element"] == "Column") & (hinge_y["Bay"] == ybay) & (
                                    hinge_y["Storey"] == st)].reset_index(drop=True).iloc[0]
                            elements["Columns"]["y"].append(et)

                        else:
                            # Columns of gravity frames
                            eleHinge = hinge_gr[(hinge_gr["Element"] == "Column") & (hinge_gr["Storey"] ==
                                                                                     st)].reset_index(drop=True).iloc[0]
                            elements["Columns"]["gravity"].append(et)

                        # Base columns
                        if st == 1:
                            base_cols.append(et)

                        self.file.write(f"\n# Column {et}, connected by {inode}, {jnode} nodes:")
                        s.hysteretic_hinges(et, inode, jnode, eleHinge, transfTag, self.flag3d, self.file)

            # Add beam elements in X direction
            for ybay in range(1, int(nbays_y + 2)):
                for xbay in range(1, int(nbays_x + 1)):
                    next_bay_x = xbay + 1
                    for st in range(1, int(nst + 1)):
                        # Element and node tags
                        et = int(f"3{xbay}{ybay}{st}")
                        inode = int(f"{xbay}{ybay}{st}")
                        jnode = int(f"{next_bay_x}{ybay}{st}")
                        # Transformation tag
                        transfTag = self.BEAM_X_TRANSF_TAG

                        if ybay == 1 or ybay == nbays_y + 1:
                            # Beams within seismic frames in X direction
                            eleHinge = hinge_x[(hinge_x["Element"] == "Beam") & (hinge_x["Bay"] == xbay) & (
                                    hinge_x["Storey"] == st)].reset_index(drop=True).iloc[0]
                            elements["Beams"]["x"].append(et)
                        else:
                            # Beams within gravity frames in X direction
                            eleHinge = hinge_gr[(hinge_gr["Element"] == "Beam") & (hinge_gr["Storey"] == st)
                                                & (hinge_gr["Direction"] == 0)].reset_index(drop=True).iloc[0]
                            elements["Beams"]["gravity_x"].append(et)

                        self.file.write(f"\n# Beam in X {et}, connected by {inode}, {jnode} nodes:")
                        s.hysteretic_hinges(et, inode, jnode, eleHinge, transfTag, self.flag3d, self.file)

            # Add beam elements in Y direction
            for xbay in range(1, int(nbays_x + 2)):
                for ybay in range(1, int(nbays_y + 1)):
                    next_bay_y = ybay + 1
                    for st in range(1, int(nst + 1)):
                        # Element and node tags
                        et = int(f"2{xbay}{ybay}{st}")
                        inode = int(f"{xbay}{ybay}{st}")
                        jnode = int(f"{xbay}{next_bay_y}{st}")

                        # Transformation tag
                        transfTag = self.BEAM_Y_TRANSF_TAG

                        if xbay == 1 or xbay == nbays_x + 1:
                            eleHinge = hinge_y[(hinge_y["Element"] == "Beam") & (hinge_y["Bay"] == ybay) & (
                                    hinge_y["Storey"] == st)].reset_index(drop=True).iloc[0]
                            elements["Beams"]["y"].append(et)
                        else:
                            eleHinge = hinge_gr[(hinge_gr["Element"] == "Beam") & (hinge_gr["Storey"] == st)
                                                & (hinge_gr["Direction"] == 1)].reset_index(drop=True).iloc[0]
                            elements["Beams"]["gravity_y"].append(et)

                        self.file.write(f"\n# Beam in Y {et}, connected by {inode}, {jnode} nodes:")
                        s.hysteretic_hinges(et, inode, jnode, eleHinge, transfTag, self.flag3d, self.file)

        else:
            # Initialize elements
            elements = {'Columns': [], 'Beams': []}

            for index, ele in self.sections.iterrows():
                if ele['Element'] == 'Beam':
                    # Beam elements
                    et = f"2{ele['Bay']}{ele['Storey']}"
                    transfTag = self.BEAM_X_TRANSF_TAG
                    elements['Beams'].append(et)

                else:
                    # Column elements
                    et = f"1{ele['Bay']}{ele['Storey']}"
                    transfTag = self.COL_TRANSF_TAG
                    elements['Columns'].append(et)
                    if ele['Storey'] == 1:
                        base_cols.append(et)

                s.hysteretic_hinges(et, None, None, ele, transfTag, self.flag3d, self.file)

        return elements, base_cols

    def rigid_diaphragm(self):
        self.file.write("\n\n# Construct a rigid diaphragm")

        nbays_x = max(self.sections["x"]["Bay"] - 1)
        nbays_y = max(self.sections["y"]["Bay"] - 1)
        # Define Rigid floor diaphragm
        master_nodes = []
        cnt = 0
        for st in range(1, self.g.nst + 1):
            hive = int(f"{int(nbays_x / 2 + 1)}{int(nbays_y / 2 + 1)}{st}")
            master_nodes.append(hive)
            # Define rigid diaphragm
            for xbay in range(int(nbays_x + 1)):
                for ybay in range(int(nbays_y + 1)):
                    bee = int(f"{1 + xbay}{1 + ybay}{st}")
                    if bee != hive:
                        op.rigidDiaphragm(3, hive, bee)
                        self.file.write(f"\nrigidDiaphragm 3 {hive} {bee};")
            cnt += 1

    def model(self):
        """
        Creates the full model
        :return: None
        """
        self.create_model()
        self.define_transformations()
        self.base_nodes, hinge_nodes = self.create_nodes()
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

        if not self.flag3d:
            # Required only for 2D modelling
            self.define_pdelta_columns(option='EqualDOF')
        else:
            self.rigid_diaphragm()

        self.define_masses()
