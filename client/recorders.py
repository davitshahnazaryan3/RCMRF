"""
Defines recorders for the model
Note: Recorders are very much generic, more options may be added as it suits the user
"""
import openseespy.opensees as op
import numpy as np


class Recorders:
    def __init__(self, geometry, elements, hingeModel='hysteretic', flag3d=False):
        """
        Initializing recorder generation
        :param geometry: Object                     Geometry object
        :param elements: dict(list(int))            Element IDs
        :param hingeModel: str                      Hinge Model
        :param flag3d: bool
        """
        self.geometry = geometry
        self.elements = elements
        self.hingeModel = hingeModel
        self.flag3d = flag3d

    def st_recorder(self, base_nodes):
        """
        Recorder for static analysis
        :param base_nodes: list(int)                Base node IDs
        :return: dict                               Dictionary including all recorders of interest
        """
        results = {"Nodes": {}}
        op.reactions()
        # Node recorders
        for n in range(len(base_nodes)):
            if self.hingeModel == 'haselton':
                node = int(f"{base_nodes[n]}0")
            else:
                node = base_nodes[n]
            results["Nodes"][node] = op.nodeReaction(node)

        # Element recorders
        results['Element'] = {}
        results['Element']['Beam'] = {}
        results['Element']['Column'] = {}

        if self.flag3d:
            for d in self.elements['Beams']:
                for beam in self.elements['Beams'][d]:
                    results["Element"]["Beam"][beam] = op.eleForce(beam)
        else:
            for beam in self.elements['Beams']:
                results['Element']['Beam'][beam] = op.eleForce(int(beam))

        if self.hingeModel == 'haselton':
            for cols in self.elements['Columns internal']:
                for ele in self.elements['Columns internal'][cols]:
                    results['Element']['Column'][ele] = op.eleForce(ele)

            for cols in self.elements['Columns external']:
                for ele in self.elements['Columns external'][cols]:
                    results['Element']['Column'][ele] = op.eleForce(ele)

        else:
            if self.flag3d:
                for d in self.elements["Columns"]:
                    for col in self.elements["Columns"][d]:
                        results["Element"]["Column"][col] = op.eleForce(col)
            else:
                for col in self.elements['Columns']:
                    results['Element']['Column'][col] = op.eleForce(int(col))

        return results

    def ma_recorder(self, num_modes, lam, path):
        """
        Records modal shapes
        :param num_modes: int                       Number of modal shapes to record
        :param lam: list                            Eigenvectors
        :param path: str                            Path to export Modal_analysis.tcl to
        :return: dict                               Dictionary containing modal shape information
        """
        if self.flag3d:

            # Get all node rags
            nodes = op.getNodeTags()

            # Initialize mass computation
            total_mass = np.array([0] * 6)

            # Compute total masses
            masses = np.zeros(self.geometry.nst)

            for node in nodes:
                indf = len(op.nodeDisp(node))
                for i in range(indf):
                    total_mass[i] += op.nodeMass(node, i + 1)

                node = str(node)
                if node[-1] != "0":
                    masses[int(node[-1]) - 1] += op.nodeMass(int(node), 1)

            # Results for each mode
            mode_data = np.zeros((num_modes, 4))
            mode_MPM = np.zeros((num_modes, 6))
            mode_L = np.zeros((num_modes, 6))

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
                L = np.zeros((6,))
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
                MPM = np.zeros((6,))
                for i in range(6):
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

            results = {"Mode1": [], "Mode2": []}

            # Initialize modal shape
            modalShape = np.zeros((self.geometry.nst, 2))
            for st in range(self.geometry.nst):
                nodetag = int(f"{self.geometry.nbays[0] + 1}{self.geometry.nbays[1] + 1}{st + 1}")

                # Mode 1 refers to X direction, and Mode 2 refers to Y direction
                results["Mode1"].append(op.nodeEigenvector(nodetag, 1, int(positions[0] + 1)))
                results["Mode2"].append(op.nodeEigenvector(nodetag, 2, int(positions[1] + 1)))
                # file.write(f"\nlappend mode1 [nodeEigenvector {nodetag} 1 {int(positions[0] + 1)}]")
                # file.write(f"\nlappend mode2 [nodeEigenvector {nodetag} 2 {int(positions[1] + 1)}]")

                # First mode shape (also for 2D model)
                modalShape[st, 0] = op.nodeEigenvector(nodetag, 1, int(positions[0] + 1))
                # Second mode shape
                modalShape[st, 1] = op.nodeEigenvector(nodetag, 2, int(positions[1] + 1))

            # Normalize the modal shapes (first two modes, most likely associated with X and Y directions unless
            # there are large torsional effects)
            modalShape = np.abs(modalShape) / np.max(np.abs(modalShape), axis=0)

            # Calculate the first mode participation factor and effective modal mass
            M = np.zeros((self.geometry.nst, self.geometry.nst))
            for st in range(self.geometry.nst):
                M[st][st] = masses[st]

            # Identity matrix
            identity = np.ones((1, self.geometry.nst))

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

            if path:
                file = open(path / "Models/modal_recorders.tcl", "w+")
                file.write("# Extracting first two modal shapes")
                # file.write("\nset mode1 {};")
                # file.write("\nset mode2 {};")

                nstart = int(f"{self.geometry.nbays[0] + 1}{self.geometry.nbays[1] + 1}{1}")
                nend = int(f"{self.geometry.nbays[0] + 1}{self.geometry.nbays[1] + 1}{self.geometry.nst}")

                file.write('\nrecorder Node -file mode1.txt -nodeRange ' +
                           f'{nstart} {nend} -dof {int(positions[0] + 1)} "eigen 1";')
                file.write('\nrecorder Node -file mode2.txt -nodeRange ' +
                           f'{nstart} {nend} -dof {int(positions[1] + 1)} "eigen 2";')
                file.close()

        else:
            results = {}
            for k in range(min(num_modes, 3)):
                results[f"Mode{k + 1}"] = []
                for st in range(self.geometry.nst):
                    if self.hingeModel == 'haselton':
                        results[f"Mode{k + 1}"].append(op.nodeEigenvector(int(f"{st + 2}{self.geometry.nbays + 1}1"),
                                                                          k + 1, 1))
                    else:
                        results[f"Mode{k + 1}"].append(op.nodeEigenvector(int(f"{self.geometry.nbays + 1}{st + 1}"),
                                                                          k + 1, 1))
            positions = None

        return results, positions
