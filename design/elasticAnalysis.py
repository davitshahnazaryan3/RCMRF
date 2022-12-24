from client.geometry import Geometry
from postprocess.opsv import *
import openseespy.opensees as op


class ElasticAnalysis:
    def __init__(self, spans_x, spans_y, heights, cross_sections, concrete_modulus, fstiff, loads=None, flag3d=False):
        """

        Parameters
        ----------
        spans_x
        spans_y
        heights
        cross_sections
        concrete_modulus
        fstiff
        loads
        flag3d: bool
            True for 3D, False for 2D
        """
        self.spans_x = spans_x
        self.spans_y = spans_y
        self.heights = heights
        self.cross_sections = cross_sections
        self.concrete_modulus = concrete_modulus
        self.fstiff = fstiff
        self.loads = loads
        self.flag3d = flag3d

        # Geometric transformation tags
        self.COL_TRANSF_TAG = 1
        self.BEAM_X_TRANSF_TAG = 2
        self.BEAM_Y_TRANSF_TAG = 3

        self.NEGLIGIBLE = 1.e-9
        self.UBIG = 1.e10

        self.base_cols = []

    def create_model(self):
        op.wipe()

        if self.flag3d:
            op.model('Basic', '-ndm', 3, '-ndf', 6)
        else:
            op.model('Basic', '-ndm', 2, '-ndf', 3)

    def define_nodes(self):
        spans_x = self.spans_x
        spans_y = self.spans_y
        heights = self.heights

        ybay = 0
        xloc = 0.
        for xbay in range(len(spans_x) + 1):
            yloc = 0.
            zloc = 0.

            # Generate the remaining nodes
            for st in range(len(heights) + 1):
                if st == 0:
                    zloc = 0.0
                else:
                    zloc += heights[st - 1]

                if self.flag3d:
                    node = int(f"{1 + xbay}{1 + ybay}{st}")
                    op.node(node, xloc, yloc, zloc)
                else:
                    node = int(f"{1 + xbay}{st}")
                    op.node(node, xloc, zloc)

                if zloc == 0:
                    if self.flag3d:
                        op.fix(node, 1, 1, 1, 1, 1, 1)
                    else:
                        op.fix(node, 1, 1, 1)

                if self.flag3d:
                    for ybay in range(len(spans_y) + 1):
                        # Base nodes
                        node = int(f"{1+xbay}{1+ybay}{st}")
                        op.node(node, xloc, yloc, zloc)

                        # Fix all base nodes
                        if zloc == 0:
                            op.fix(node, 1, 1, 1, 1, 1, 1)

                        # increment y coordinate
                        if ybay != len(spans_y):
                            yloc += spans_y[ybay]

            # increment x coordinate
            if xbay != len(spans_x):
                xloc += spans_x[xbay]

    def define_geometric_transformations(self, name="Linear"):
        if self.flag3d:
            op.geomTransf(name, self.COL_TRANSF_TAG, 0., -1., 0., '-jntOffset', 0.0, 0.0, 0.3, 0.0, 0.0, -0.3)
            op.geomTransf(name, self.BEAM_X_TRANSF_TAG, 0., -1., 0., '-jntOffset', 0.2, 0.0, 0.0, -0.2, 0.0, 0.0)
            op.geomTransf(name, self.BEAM_Y_TRANSF_TAG, 1., 0., 0., '-jntOffset', 0.0, 0.2, 0.0, 0.0, -0.2, 0.0)
        else:
            op.geomTransf(name, self.COL_TRANSF_TAG, 0., -1., 0., '-jntOffset', 0.0, 0.3, 0.0, -0.3)
            op.geomTransf(name, self.BEAM_X_TRANSF_TAG, 0., -1., 0., '-jntOffset', 0.2, 0.0, -0.2, 0.0)

    def create_elements(self):

        cs_x = self.cross_sections["x_seismic"]

        spans_x = self.spans_x
        heights = self.heights

        if not self.flag3d:
            columns = []
            beams = []

            # Columns
            for xbay in range(1, len(spans_x) + 2):
                for st in range(1, len(heights) + 1):
                    # previous storey level
                    previous_st = st - 1

                    et = int(f"1{xbay}{st}")
                    columns.append(et)

                    if xbay == 1 or xbay == len(spans_x) + 1:
                        b = h = cs_x[f"he{st}"]
                    else:
                        b = h = cs_x[f"hi{st}"]

                    area = b * h
                    iz = b * h**3 / 12

                    # End nodes of column
                    inode = int(f"{xbay}{previous_st}")
                    jnode = int(f"{xbay}{st}")

                    op.element("elasticBeamColumn", et, inode, jnode, area, self.concrete_modulus * self.fstiff,
                               iz, self.COL_TRANSF_TAG)

            # Beams
            for st in range(1, len(heights) + 1):
                for xbay in range(1, len(spans_x) + 1):
                    next_bay_x = xbay + 1

                    # element and node tags
                    et = int(f"2{xbay}{st}")
                    beams.append(et)

                    b = cs_x[f"b{st}"]
                    h = cs_x[f"h{st}"]

                    area = b * h
                    iz = b * h**3 / 12

                    # End nodes of column
                    inode = int(f"{xbay}{st}")
                    jnode = int(f"{next_bay_x}{st}")

                    op.element("elasticBeamColumn", et, inode, jnode, area, self.concrete_modulus * self.fstiff,
                               iz, self.BEAM_X_TRANSF_TAG)

            return beams, columns

        # For storing
        columns = {"x": [], "y": [], "gravity": []}
        beams = {"x": [], "y": [], "gravity_x": [], "gravity_y": []}

        cs_y = self.cross_sections["y_seismic"]
        cs_gr = self.cross_sections["gravity"]

        spans_y = self.spans_y
        Gc = self.concrete_modulus / 2.0 / (1 + 0.2)

        # Add column elements
        for xbay in range(1, len(spans_x) + 2):
            for st in range(1, len(heights) + 1):
                # previous storey level
                previous_st = st - 1

                for ybay in range(1, len(spans_y) + 2):
                    # Column element tag
                    et = int(f"1{xbay}{ybay}{st}")

                    b, h = self.get_column_cross_section_dimensions(xbay, ybay, st, cs_x, cs_y, cs_gr)
                    area = b * h
                    iz = b * h ** 3 / 12
                    iy = h * b ** 3 / 12

                    # Torsional moment of inertia
                    if h >= b:
                        J = b * h ** 3 * (16 / 3 - 3.36 * h / b * (1 - 1 / 12 * (h / b) ** 4))
                    else:
                        J = h * b ** 3 * (16 / 3 - 3.36 * b / h * (1 - 1 / 12 * (b / h) ** 4))

                    # End nodes of column
                    inode = int(f"{xbay}{ybay}{previous_st}")
                    jnode = int(f"{xbay}{ybay}{st}")

                    # Create elasticBeamColumn
                    op.element("elasticBeamColumn", et, inode, jnode, area, self.concrete_modulus * self.fstiff, Gc, J,
                               iy, iz, self.COL_TRANSF_TAG)
                    # For recorders
                    if ybay == 1 or ybay == len(spans_y) + 1:
                        columns["x"].append(et)
                    elif (xbay == 1 or xbay == len(spans_x) + 1) and (1 < ybay < len(spans_y) + 1):
                        columns["y"].append(et)
                    else:
                        columns["gravity"].append(et)

                    # Base columns (for recorders, because for some reason base node recorders did not record)
                    if st == 1:
                        self.base_cols.append(et)

        # Add beam elements in X direction
        for ybay in range(1, len(spans_y) + 2):
            for st in range(1, len(heights) + 1):
                for xbay in range(1, len(spans_x) + 1):
                    next_bay_x = xbay + 1

                    # element and node tags
                    et = int(f"3{xbay}{ybay}{st}")
                    inode = int(f"{xbay}{ybay}{st}")
                    jnode = int(f"{next_bay_x}{ybay}{st}")

                    b, h = self.get_beam_cross_section_dimensions("x", ybay, len(spans_y), st, cs_x, cs_gr)
                    area = b * h
                    iz = b * h ** 3 / 12
                    iy = h * b ** 3 / 12

                    # Torsional moment of inertia
                    if h >= b:
                        J = b * h ** 3 * (16 / 3 - 3.36 * h / b * (1 - 1 / 12 * (h / b) ** 4))
                    else:
                        J = h * b ** 3 * (16 / 3 - 3.36 * b / h * (1 - 1 / 12 * (b / h) ** 4))

                    # Create elasticBeamColumn
                    op.element("elasticBeamColumn", et, inode, jnode, area, self.concrete_modulus * self.fstiff, Gc, J,
                               iy, iz, self.BEAM_X_TRANSF_TAG)

                    # For recorders
                    if ybay == 1 or ybay == len(spans_y) + 1:
                        beams["x"].append(et)
                    else:
                        beams["gravity_x"].append(et)

        # Add beam elements in Y direction
        for xbay in range(1, len(spans_x) + 2):
            for ybay in range(1, len(spans_y) + 1):
                next_bay_y = ybay + 1
                for st in range(1, len(heights) + 1):
                    # element and node tags
                    et = int(f"2{xbay}{ybay}{st}")
                    inode = int(f"{xbay}{ybay}{st}")
                    jnode = int(f"{xbay}{next_bay_y}{st}")

                    b, h = self.get_beam_cross_section_dimensions("y", xbay, len(spans_x), st, cs_y, cs_gr)
                    area = b * h
                    iz = b * h ** 3 / 12
                    iy = h * b ** 3 / 12

                    # Torsional moment of inertia
                    if h >= b:
                        J = b * h ** 3 * (16 / 3 - 3.36 * h / b * (1 - 1 / 12 * (h / b) ** 4))
                    else:
                        J = h * b ** 3 * (16 / 3 - 3.36 * b / h * (1 - 1 / 12 * (b / h) ** 4))

                    # Create elasticBeamColumn
                    op.element("elasticBeamColumn", et, inode, jnode, area, self.concrete_modulus * self.fstiff, Gc, J,
                               iy, iz, self.BEAM_Y_TRANSF_TAG)

                    # For recorders
                    if xbay == 1 or xbay == len(spans_x) + 1:
                        beams["y"].append(et)
                    else:
                        beams["gravity_y"].append(et)

        return beams, columns

    @staticmethod
    def get_beam_cross_section_dimensions(direction, bay, nbays, st, cs, cs_gr):
        if bay == 1 or bay == nbays + 1:
            b = cs[f"b{st}"]
            h = cs[f"h{st}"]
        else:
            b = cs_gr[f"b{direction}{st}"]
            h = cs_gr[f"h{direction}{st}"]

        return b, h

    def get_column_cross_section_dimensions(self, xbay, ybay, st, cs_x, cs_y, cs_gr):
        if ybay == 1 or ybay == len(self.spans_y) + 1:
            if xbay == 1 or xbay == len(self.spans_x) + 1:
                # External columns
                b = h = cs_x[f"he{st}"]
            else:
                # Interior columns
                b = h = cs_x[f"hi{st}"]

        # Columns of seismic frame along y direction
        elif (xbay == 1 or xbay == len(self.spans_x) + 1) and (1 < ybay < len(self.spans_y) + 1):
            # External columns are already created
            # Internal columns
            b = h = cs_y[f"hi{st}"]

        # Columns of gravity frames
        else:
            # Only internal columns
            b = h = cs_gr[f"hi{st}"]

        return b, h

    def define_gravity_loads(self, beams, grav_loads=None):
        if not self.flag3d:
            # assuming that external spans are equal in length, and we have perimeter seismic frames
            # additionally: the load is assumed to be transferred to X beams only (conservative)
            loads = np.array(self.loads["loads"]) * self.spans_y[0] / 2
            grav_loads = loads

        Ew = {}

        op.timeSeries('Constant', 1)
        op.pattern('Plain', 1, 1)

        if grav_loads is not None and None not in grav_loads:
            if self.flag3d:
                # Seismic frames
                for ele in beams["x"]:
                    st = int(str(ele)[-1]) - 1
                    op.eleLoad('-ele', ele, '-type', '-beamUniform', -abs(grav_loads["x"][st]), self.NEGLIGIBLE)
                    Ew[ele] = ['-beamUniform', -abs(grav_loads["x"][st]), self.NEGLIGIBLE, self.NEGLIGIBLE]
                for ele in beams["y"]:
                    st = int(str(ele)[-1]) - 1
                    op.eleLoad('-ele', ele, '-type', '-beamUniform', -abs(grav_loads["y"][st]), self.NEGLIGIBLE)
                    Ew[ele] = ['-beamUniform', -abs(grav_loads["y"][st]), self.NEGLIGIBLE, self.NEGLIGIBLE]
                # Gravity frames
                for ele in beams["gravity_x"]:
                    st = int(str(ele)[-1]) - 1
                    op.eleLoad('-ele', ele, '-type', '-beamUniform', -2 * abs(grav_loads["x"][st]), self.NEGLIGIBLE)
                    Ew[ele] = ['-beamUniform', -2*abs(grav_loads["x"][st]), self.NEGLIGIBLE, self.NEGLIGIBLE]
                for ele in beams["gravity_y"]:
                    st = int(str(ele)[-1]) - 1
                    op.eleLoad('-ele', ele, '-type', '-beamUniform', -2 * abs(grav_loads["y"][st]), self.NEGLIGIBLE)
                    Ew[ele] = ['-beamUniform', -2*abs(grav_loads["y"][st]), self.NEGLIGIBLE, self.NEGLIGIBLE]
            else:
                for ele in beams:
                    storey = int(str(ele)[-1]) - 1
                    op.eleLoad('-ele', ele, '-type', '-beamUniform', -abs(grav_loads[storey]))
                    Ew[ele] = ['-beamUniform', -abs(grav_loads[storey]), 0.0]

        else:
            Ew = self.compute_gravity_loads(beams, Ew)

        return Ew

    def compute_gravity_loads(self, beams, Ew=None, distributed=True):
        if Ew is None:
            Ew = {}

        # get distributed loads and span dimensions
        q_roof = self.loads['loads'][1]
        q_floor = self.loads['loads'][0]
        spans_x = self.spans_x
        spans_y = self.spans_y

        for d in beams:
            for beam in beams[d]:
                st = int(str(beam)[-1])
                xbay = int(str(beam)[1])
                ybay = int(str(beam)[2])

                # Distributed loads (kN/m2)
                if st != len(self.heights):
                    # General storey
                    q = q_floor
                else:
                    # Roof level
                    q = q_roof

                # Distributing the loads
                if d == "x" or d == "gravity_x":
                    # Beams along X direction
                    # Load over a beam
                    control_length = spans_y[ybay - 1] if ybay < len(spans_y) + 1 else spans_y[ybay - 2]
                    Ew = self.get_load(beam, spans_x, control_length, xbay, q, distributed, Ew=Ew)

                    # Additional load for interior beams
                    if 1 < ybay < len(spans_y) + 1:
                        control_length = spans_y[ybay - 2]
                        Ew = self.get_load(beam, spans_x, control_length, xbay, q, distributed, Ew=Ew)

                else:
                    # Beams along Y direction
                    # Load over a beam
                    control_length = spans_x[xbay - 1] if xbay < len(spans_x) + 1 else spans_x[xbay - 2]
                    Ew = self.get_load(beam, spans_y, control_length, ybay, q, distributed, 1, Ew=Ew)

                    # Additional load for interior beams
                    if 1 < xbay < len(spans_x) + 1:
                        control_length = spans_x[xbay - 2]
                        Ew = self.get_load(beam, spans_y, control_length, ybay, q, distributed, 1, Ew=Ew)
        return Ew

    def get_load(self, beam, spans, control_length, bay, q, distributed, direction=0, Ew=None):
        if spans[bay - 1] <= control_length:
            # Triangular rule
            load = q * spans[bay - 1] ** 2 / 4 / spans[bay - 1]
        else:
            # Trapezoidal rule
            load = 1 / 4 * q * control_length * (2 * spans[bay - 1] - control_length) / \
                   spans[bay - 1]
        load = round(load, 2)

        # apply the load
        if distributed:
            # If the load is applied uniformly
            op.eleLoad('-ele', beam, '-type', '-beamUniform', load, self.NEGLIGIBLE)
            if beam in Ew:
                Ew[beam][1] += load

            else:
                Ew[beam] = ['-beamUniform', -load, self.NEGLIGIBLE, self.NEGLIGIBLE]

        else:
            # End nodes connecting the beam
            if direction == 0:
                inode = beam - 3000
                jnode = beam - 3000 + 100
            else:
                inode = beam - 2000
                jnode = beam - 2000 + 10

            # If the load is applied as point load at the end nodes of the beam
            op.load(inode, self.NEGLIGIBLE, self.NEGLIGIBLE, -load * spans[bay - 1] / 2,
                    self.NEGLIGIBLE, self.NEGLIGIBLE, self.NEGLIGIBLE)
            op.load(jnode, self.NEGLIGIBLE, self.NEGLIGIBLE, -load * spans[bay - 1] / 2,
                    self.NEGLIGIBLE, self.NEGLIGIBLE, self.NEGLIGIBLE)

        return Ew

    def define_lateral_loads(self, action, direction=0):
        spans_x = self.spans_x
        spans_y = self.spans_y
        heights = self.heights

        op.timeSeries("Constant", 2)
        op.pattern("Plain", 2, 2)

        if not self.flag3d:
            n_nodes_x = len(spans_x) + 1
            n_nodes_y = len(spans_y) + 1
        else:
            n_nodes_x = n_nodes_y = (len(spans_y) + 1) * (len(spans_x) + 1)

        for st in range(1, len(heights) + 1):
            if direction == 0:
                # Along x direction
                for bay in range(1, len(spans_x) + 2):
                    self.apply_lateral(action, len(spans_y), n_nodes_x, bay, st)
            else:
                # Along y direction
                for bay in range(1, len(spans_y) + 2):
                    self.apply_lateral(action, len(spans_x), n_nodes_y, bay, st, direction)

    def apply_lateral(self, lateral_action, nbay, nodes, bay, st, direction=0):

        if direction == 0:
            x = lateral_action[st - 1] / nodes
            y = self.NEGLIGIBLE

            # element IDs
            ele1 = int(f"{bay}1{st}")
            ele2 = int(f"{bay}{nbay + 1}{st}")
        else:
            x = self.NEGLIGIBLE
            y = lateral_action[st - 1] / nodes

            # element IDs
            ele1 = int(f"1{bay}{st}")
            ele2 = int(f"{nbay + 1}{bay}{st}")

        if not self.flag3d:
            # 2D model
            op.load(int(f"{bay}{st}"), y, self.NEGLIGIBLE, self.NEGLIGIBLE)
            return

        op.load(ele1, x, y, self.NEGLIGIBLE, self.NEGLIGIBLE, self.NEGLIGIBLE, self.NEGLIGIBLE)
        op.load(ele2, x, y, self.NEGLIGIBLE, self.NEGLIGIBLE, self.NEGLIGIBLE, self.NEGLIGIBLE)

        for ybay in range(2, int(nbay + 1)):
            if direction == 0:
                ele = int(f"{bay}{ybay}{st}")
            else:
                ele = int(f"{ybay}{bay}{st}")

            op.load(ele, x, y, self.NEGLIGIBLE, self.NEGLIGIBLE, self.NEGLIGIBLE, self.NEGLIGIBLE)

    def run_static_analysis(self):
        if self.flag3d:
            dgravity = 1.0 / 1
            op.integrator("LoadControl", dgravity)
            op.numberer("RCM")
            op.system("UmfPack")
            op.constraints("Penalty", 1.0e15, 1.0e15)
            op.test("EnergyIncr", 1.0e-8, 10)
        else:
            op.integrator("LoadControl", 0.1)
            op.numberer("Plain")
            op.system("BandGeneral")
            op.constraints("Plain")
            op.test("NormDispIncr", 1.0e-8, 6)
        op.algorithm("Newton")
        op.analysis("Static")
        op.analyze(1)
        op.loadConst("-time", 0.0)

    def get_section_force_diagrams(self, Ew):
        if self.flag3d:
            forces = section_force_diagram_3d(Ew)
        else:
            forces = section_force_diagram_2d(Ew)
        return forces

    def run_elastic_analysis(self, gravity=True, lateral=False, direction=0, lat_action=None, gravity_loads=None):
        self.create_model()
        self.define_nodes()
        self.define_geometric_transformations()
        beams, columns = self.create_elements()

        if gravity:
            Ew = self.define_gravity_loads(beams, gravity_loads)

        if lateral:
            self.define_lateral_loads(lat_action, direction=direction)

        self.run_static_analysis()

        # Define recorders
        demands = self.record(len(self.spans_x), len(self.spans_y), direction=direction)

        diagrams = None
        if gravity:
            temp = self.get_section_force_diagrams(Ew)
            diagrams = temp.copy()

        op.wipe()

        return demands, diagrams

    def record(self, nbays_x, nbays_y, direction=0):
        nst = len(self.heights)

        # Indices for recorders
        if self.flag3d:
            if direction == 0:
                # Along X axis
                midx = [11, 5]
                vidx_c = [1, 7]
                vidx_b = [3, 9]
                nidx_c = [3, 9]
                nidx_b = [1, 7]
            else:
                # Along Y axis
                midx = [10, 4]
                vidx_c = [2, 8]
                vidx_b = [3, 9]
                nidx_c = [3, 9]
                nidx_b = [2, 8]
        else:
            midx = [3, 6]
            vidx_c = [1, 4]
            vidx_b = [2, 5]
            nidx_c = [2, 5]
            nidx_b = [1, 4]

        # Define recorders
        # Seismic frame element recorders initialization
        # Only bx_seismic and cx_seismic are of interest for 2D modelling
        bx_seismic = np.zeros((nst, nbays_x))
        by_seismic = np.zeros((nst, nbays_y))
        cx_seismic = np.zeros((nst, nbays_x + 1))
        cy_seismic = np.zeros((nst, nbays_y + 1))

        # For gravity frame elements only the max demands will be used for uniform design
        bx_gravity = np.zeros((nst, nbays_x, nbays_y - 1))
        by_gravity = np.zeros((nst, nbays_x - 1, nbays_y))
        c_gravity = np.zeros((nst, nbays_x - 1, nbays_y - 1))

        # Global dictionary variable to store all demands
        # Only x_seismic is relevant for 2D modelling
        results = {"x_seismic": {"Beams": {"M": {"Pos": bx_seismic.copy(), "Neg": bx_seismic.copy()},
                                           "N": bx_seismic.copy(), "V": bx_seismic.copy()},
                                 "Columns": {"M": cx_seismic.copy(), "N": cx_seismic.copy(), "V": cx_seismic.copy()}},
                   "y_seismic": {"Beams": {"M": {"Pos": by_seismic.copy(), "Neg": by_seismic.copy()},
                                           "N": by_seismic.copy(), "V": by_seismic.copy()},
                                 "Columns": {"M": cy_seismic.copy(), "N": cy_seismic.copy(), "V": cy_seismic.copy()}},
                   "gravity": {"Beams_x": {"M": {"Pos": bx_gravity.copy(), "Neg": bx_gravity.copy()},
                                           "N": bx_gravity.copy(), "V": bx_gravity.copy()},
                               "Beams_y": {"M": {"Pos": by_gravity.copy(), "Neg": by_gravity.copy()},
                                           "N": by_gravity.copy(), "V": by_gravity.copy()},
                               "Columns": {"M": c_gravity.copy(), "N": c_gravity.copy(), "V": c_gravity.copy()}}}

        # Beams, counting: bottom to top, left to right
        # Beams of seismic frame along X direction
        # Beam iNode [Fx, Fy, Fz, Mx, My, Mz]; jNode [Fx, Fy, Fz, Myx, My, Mz]
        # iNode Negative My means upper demand; jNode Positive My means upper demand
        for bay in range(nbays_x):
            for st in range(nst):
                if self.flag3d:
                    et = int(f"3{bay+1}1{st+1}")
                else:
                    et = int(f"2{bay+1}{st+1}")

                results = self.record_results(results, "x_seismic", "Beams", st, bay, et, midx, nidx_b, vidx_b)

        # Columns
        # Columns of seismic frame along x direction
        # Columns [Vx, Vy, N, Mx, My, Mz]; jNode [Vx, Vy, N, Mx, My, Mz]
        # Columns for X direction demand estimations [V, 2, N, 3, M, 6]
        # Columns for Y direction demand estimations [0, V, N, M, 5, 6]
        for bay in range(nbays_x + 1):
            for st in range(nst):
                if self.flag3d:
                    et = int(f"1{bay + 1}1{st + 1}")
                else:
                    et = int(f"1{bay+1}{st+1}")

                results = self.record_results(results, "x_seismic", "Columns", st, bay, et, midx, nidx_c, vidx_c)

        if not self.flag3d:
            # return if 2D modelling was selected, otherwise continue recording
            return results["x_seismic"]

        # Remaining Beams
        # Beams of seismic frame along Y direction
        for bay in range(nbays_y):
            for st in range(nst):
                et = int(f"21{bay+1}{st+1}")
                results = self.record_results(results, "y_seismic", "Beams", st, bay, et, midx, nidx_b, vidx_b)

        # Beams of gravity frames along X direction
        for ybay in range(nbays_y - 1):
            for xbay in range(nbays_x):
                for st in range(nst):
                    et = int(f"3{xbay+1}{ybay+2}{st+1}")
                    results = self.record_results(results, "gravity", "Beams_x", st, xbay, et, midx, nidx_b, vidx_b,
                                                  ybay)

        # Beams of gravity frames along Y direction
        for xbay in range(nbays_x - 1):
            for ybay in range(nbays_y):
                for st in range(nst):
                    et = int(f"2{xbay + 2}{ybay + 1}{st + 1}")
                    results = self.record_results(results, "gravity", "Beams_y", st, xbay, et, midx, nidx_b, vidx_b,
                                                  ybay)

        # Remaining Columns
        # Columns of seismic frame along y direction
        for bay in range(nbays_y + 1):
            for st in range(nst):
                et = int(f"11{bay+1}{st+1}")
                results = self.record_results(results, "y_seismic", "Columns", st, bay, et, midx, nidx_c, vidx_c)

        # Columns of gravity frames
        for xbay in range(nbays_x - 1):
            for ybay in range(nbays_y - 1):
                for st in range(nst):
                    et = int(f"1{xbay+2}{ybay+2}{st+1}")
                    results = self.record_results(results, "gravity", "Columns", st, xbay, et, midx, nidx_c, vidx_c,
                                                  ybay)

        return results

    @staticmethod
    def record_results(results, frame, element, st, bay, et, moment, axial, shear, ybay=None):

        if ybay:
            # Interior frames
            if element == "Columns":
                results[frame][element]["M"][st][bay][ybay] = max(abs(op.eleForce(et, moment[1])),
                                                                  abs(op.eleForce(et, moment[0])))
            else:
                results[frame][element]["M"]["Pos"][st][bay][ybay] = abs(op.eleForce(et, moment[1]))
                results[frame][element]["M"]["Neg"][st][bay][ybay] = abs(op.eleForce(et, moment[0]))
            results[frame][element]["N"][st][bay][ybay] = max(op.eleForce(et, axial[0]),
                                                              op.eleForce(et, axial[1]), key=abs)
            results[frame][element]["V"][st][bay][ybay] = max(abs(op.eleForce(et, shear[0])),
                                                              abs(op.eleForce(et, shear[1])))

        else:
            if element == "Columns":
                results[frame][element]["M"][st][bay] = max(abs(op.eleForce(et, moment[1])),
                                                            abs(op.eleForce(et, moment[0])))

            else:
                results[frame][element]["M"]["Pos"][st][bay] = abs(op.eleForce(et, moment[1]))
                results[frame][element]["M"]["Neg"][st][bay] = abs(op.eleForce(et, moment[0]))
            results[frame][element]["N"][st][bay] = max(op.eleForce(et, axial[0]),
                                                        op.eleForce(et, axial[1]), key=abs)
            results[frame][element]["V"][st][bay] = max(abs(op.eleForce(et, shear[0])),
                                                        abs(op.eleForce(et, shear[1])))

        return results

    def define_masses(self, loads):
        """
        Define masses. Mass should be the total mass of the building, which is then divided by the number of seismic
        frames.
        :return: None
        """
        heights = self.heights
        spans_x = self.spans_x
        spans_y = self.spans_y

        nst = len(heights)
        nbays_x = len(spans_x)
        nbays_y = len(spans_y)

        masses = [0] * nst
        previous_x = 0
        for st in range(1, nst + 1):
            for xbay in range(1, nbays_x + 2):
                if self.flag3d:
                    for ybay in range(1, nbays_y + 2):
                        node = int(f"{xbay}{ybay}{st}")

                        # Exterior nodes
                        if xbay == 1:
                            if ybay == 1:
                                # Exterior node
                                area = spans_x[xbay - 1] * spans_y[ybay - 1] / 4
                            elif ybay == nbays_y + 1:
                                # Exterior node
                                area = spans_x[xbay - 1] * spans_y[ybay - 2] / 4
                            else:
                                # Side node
                                area = spans_x[xbay - 1] * (spans_y[ybay - 2] + spans_y[ybay - 1]) / 4

                        elif xbay == nbays_x + 1:
                            if ybay == 1:
                                # Exterior node
                                area = spans_x[xbay - 2] * spans_y[ybay - 1] / 4
                            elif ybay == nbays_y + 1:
                                # Exterior node
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
                                # Interior node
                                area = (spans_x[xbay - 2] + spans_x[xbay - 1]) * (spans_y[ybay - 2] + spans_y[ybay - 1]) / 4

                        # Mass based on tributary area
                        mass = area * loads[st-1] / 9.81
                        masses[st-1] += mass
                        op.mass(node, mass, mass, mass, self.NEGLIGIBLE, self.NEGLIGIBLE, self.NEGLIGIBLE)

                else:
                    # Assumption: perimeter seismic frames, mass is divided by 2 for each seismic frame
                    tributary_length = sum(spans_y) / 2

                    if xbay == nbays_x + 1:
                        x_length = 0
                    else:
                        x_length = spans_x[xbay-1]
                    x_total = (x_length + previous_x) / 2

                    if xbay == 1 or xbay == nbays_x + 1:
                        # corner nodes
                        mass = loads[st - 1] * tributary_length * x_total / 9.81
                    else:
                        mass = loads[st - 1] * tributary_length * x_total / 9.81

                    previous_x = x_length
                    masses[st - 1] += mass
                    op.mass(int(f"{xbay}{st}"), mass, self.NEGLIGIBLE, mass, self.NEGLIGIBLE, self.NEGLIGIBLE,
                            self.NEGLIGIBLE)

        return masses

    def run_modal_analysis(self, num_modes, masses):
        """
        Runs modal analysis
        :param num_modes: DataFrame                 Design solution, cross-section dimensions
        :return: list                               Modal periods
        """
        nbays = len(self.spans_x)
        nst = len(self.heights)

        # Get all node tags
        nodes = op.getNodeTags()

        # Check problem size (2D or 3D)
        ndm = len(op.nodeCoord(nodes[0]))

        # Initialize computation of total masses
        if ndm == 3:
            # 3D building
            ndf_max = 6
            total_mass = np.array([0] * 6)
        else:
            # 2D frame
            ndf_max = 3
            total_mass = np.array([0] * 3)

        # Get the total masses
        for node in nodes:
            indf = len(op.nodeDisp(node))
            for i in range(indf):
                total_mass[i] += op.nodeMass(node, i + 1)

        # Compute the eigenvectors (solver)
        lam = None
        try:
            lam = op.eigen(num_modes)
        except:
            print("[EXCEPTION] Eigensolver failed, trying genBandArpack...")
            try:
                lam = op.eigen('-genBandArpack', num_modes)
            except:
                print("[EXCEPTION] Eigensolver failed, trying fullGenLapack...")
                try:
                    lam = op.eigen('-fullGenLapack', num_modes)
                except:
                    print("[EXCEPTION] Eigensolver failed, trying symmBandLapack...")
                    try:
                        lam = op.eigen('-symmBandLapack', num_modes)
                    except:
                        print("[EXCEPTION] Eigensolver failed.")

        # Record stuff
        op.record()

        # Results for each mode
        mode_data = np.zeros((num_modes, 4))
        mode_MPM = np.zeros((num_modes, ndf_max))
        mode_L = np.zeros((num_modes, ndf_max))

        # Extract eigenvalues to appropriate arrays
        omega = []
        freq = []
        period = []
        for m in range(num_modes):
            omega.append(np.sqrt(lam[m]))
            freq.append(np.sqrt(lam[m]) / 2 / np.pi)
            period.append(2 * np.pi / np.sqrt(lam[m]))
            mode_data[m, :] = np.array([lam[m], omega[m], freq[m], period[m]])

            # Compute L and gm
            L = np.zeros((ndf_max,))
            gm = 0
            for node in nodes:
                V = op.nodeEigenvector(node, m + 1)
                indf = len(op.nodeDisp(node))
                for i in range(indf):
                    Mi = op.nodeMass(node, i + 1)
                    Vi = V[i]
                    Li = Mi * Vi
                    gm += Vi ** 2 * Mi
                    L[i] += Li
            mode_L[m, :] = L

            # Compute MPM
            MPM = np.zeros((ndf_max,))
            for i in range(ndf_max):
                Li = L[i]
                TMi = total_mass[i]
                MPMi = Li ** 2
                if gm > 0.0:
                    MPMi = MPMi / gm
                if TMi > 0.0:
                    MPMi = MPMi / TMi * 100.0
                MPM[i] = MPMi
            mode_MPM[m, :] = MPM

        # Get modal positions based on mass participation
        positions = np.argmax(mode_MPM, axis=1)
        # Take the first two, as for symmetric structures higher modes are not so important
        positions = positions[:2]

        # Calculate the first modal shape
        modalShape = np.zeros((nst, 2))
        for st in range(nst):
            if self.flag3d:
                nodetag = int(f"{nbays + 1}1{st + 1}")
            else:
                nodetag = int(f"{nbays+1}{st+1}")
                positions[0] = 0

            # First mode shape (also for 2D model)
            modalShape[st, 0] = op.nodeEigenvector(nodetag, 1, int(positions[0] + 1))
            # Second mode shape
            modalShape[st, 1] = op.nodeEigenvector(nodetag, 2, int(positions[1] + 1))

        # Normalize the modal shapes (first two modes, most likely associated with X and Y directions unless there are
        # large torsional effects)
        modalShape = np.abs(modalShape) / np.max(np.abs(modalShape), axis=0)

        # Calculate the first mode participation factor and effective modal mass
        M = np.zeros((nst, nst))
        for st in range(nst):
            M[st][st] = masses[st]

        # Identity matrix
        identity = np.ones((1, nst))

        gamma = np.zeros(2)
        mstar = np.zeros(2)
        for i in range(2):
            # Modal participation factor
            gamma[i] = (modalShape[:, i].transpose().dot(M)).dot(identity.transpose()) / \
                       (modalShape[:, i].transpose().dot(M)).dot(modalShape[:, i])

            # Modal mass
            mstar[i] = (modalShape[:, i].transpose().dot(M)).dot(identity.transpose())

        # Modify indices of modal properties as follows:
        # index 0 = direction x
        # index 1 = direction y
        period = np.array([period[i] for i in range(len(positions))])
        gamma = np.array([gamma[i] for i in range(len(positions))])
        mstar = np.array([mstar[i] for i in range(len(positions))])

        # Wipe analysis
        op.wipe()

        return period, modalShape, gamma, mstar


if __name__ == "__main__":
    from pathlib import Path
    import pandas as pd

    spans_x = [6., 5.]
    spans_y = [6.]
    heights = [3.5, 3.]
    gravity = True
    lateral = True
    Ew = None
    direction = 0

    directory = Path.cwd().parents[0]
    csx = directory / "sample/AnconaSaAvg/solution_cache_space_x2.csv"
    csy = directory / "sample/AnconaSaAvg/solution_cache_space_y2.csv"
    csg = directory / "sample/AnconaSaAvg/solution_cache_space_gr2.csv"
    csx = pd.read_csv(csx, index_col=0).iloc[0]
    csy = pd.read_csv(csy, index_col=0).iloc[0]
    csg = pd.read_csv(csg, index_col=0).iloc[0]
    cs = {"x_seismic": csx, "y_seismic": csy, "gravity": csg}

    lateral_action = [100, 200]
    gravity_loads = {"x": [20, 15],
                     "y": [15, 10]}

    ea = ElasticAnalysis(spans_x, spans_y, heights, cs, 31000., 1.)
    demands, diagrams = ea.run_elastic_analysis(gravity, lateral, direction, lateral_action, gravity_loads)

    if gravity:
        print("Element", 3111)
        print(np.amax(diagrams[3111], axis=0))
        print(np.amin(diagrams[3111], axis=0))

        # for ele in diagrams:
        #     print("Element", ele)
        #     print(np.amax(diagrams[ele], axis=0))
        #     print(np.amin(diagrams[ele], axis=0))
