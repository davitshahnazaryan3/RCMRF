"""
Defines recorders for the model
Note: Recorders are very much generic, more options may be added as it suits the user
"""
import openseespy.opensees as op


class Recorders:
    def __init__(self, nbays, nst, elements, hingeModel='haselton'):
        """
        Initializing recorder generation
        :param nbays: int                           Number of bays
        :param nst: int                             Number of stories
        :param elements: dict(list(int))            Element IDs
        :param hingeModel: str                      Hinge Model
        """
        self.nbays = nbays
        self.nst = nst
        self.elements = elements
        self.hingeModel = hingeModel

    def st_recorder(self, base_nodes):
        """
        Recorder for static analysis
        :param base_nodes: list(int)                Base node IDs
        :return: dict                               Dictionary including all recorders of interest
        """
        results = {}
        op.reactions()
        # Node recorders
        for n in range(len(base_nodes)):
            if self.hingeModel == 'haselton':
                node = int(f"{base_nodes[n]}0")
            else:
                node = base_nodes[n]
            results[node] = op.nodeReaction(node)

        # Element recordders
        results['Element'] = {}
        results['Element']['Beam'] = {}
        results['Element']['Column'] = {}
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
            for col in self.elements['Columns']:
                results['Element']['Column'][col] = op.eleForce(int(col))

        return results

    def ma_recorder(self, num_modes):
        """
        Records modal shapes
        :param num_modes: int                       Number of modal shapes to record
        :return: dict                               Dictionary containing modal shape information
        """
        results = {}
        for k in range(min(num_modes, 3)):
            results[f"Mode{k + 1}"] = []
            for st in range(self.nst):
                if self.hingeModel == 'haselton':
                    results[f"Mode{k + 1}"].append(op.nodeEigenvector(int(f"{st + 2}{self.nbays + 1}1"), k + 1, 1))
                else:
                    results[f"Mode{k + 1}"].append(op.nodeEigenvector(int(f"{self.nbays + 1}{st + 1}"), k + 1, 1))

        return results
