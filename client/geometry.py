"""
Creates geometry of the model
"""
import pandas as pd
import numpy as np


class Geometry:
    def __init__(self, sections, hingeModel, flag3d=False, direction=0):
        """
        initializes geometry creation
        :param sections: DataFrame              Section properties
        :param hingeModel: str                  Hinge model type (Haselton (4 nodes) or Hysteretic (2 nodes))
        :param flag3d: bool
        :param direction: int
        """
        self.sections = sections
        self.hingeModel = hingeModel.lower()
        self.flag3d = flag3d
        self.direction = direction

        self.nst = None
        self.nbays = []
        self.heights = None
        self.widths = []
        self.beams = []
        self.columns = []
        self.bnode = []
        self.tnode = []

        if self.flag3d:
            for i in range(2):
                if i == 0:
                    selection = self.sections["x"]
                else:
                    selection = self.sections["y"]

                nst, nbays, heights, widths, beams, columns, bnode, tnode = self._generate_data(selection, i)

                self.nbays.append(nbays)
                self.widths.append(widths)
                self.beams.append(beams)
                self.columns.append(columns)
                self.bnode.append(bnode)
                self.tnode.append(tnode)
                self.nst = nst
                self.heights = heights

        else:
            self.nst, self.nbays, self.heights, self.widths, self.beams, self.columns, self.bnode, self.tnode = \
                self._generate_data(self.sections, direction=self.direction)

    def _generate_data(self, sections, direction):
        # Get the number of storeys and bays in the direction of seismic action
        if self.hingeModel == 'haselton':
            nst = max(sections['Storey'])
            nbays = int(len(sections[(sections['Element'] == 'Column')]) / nst - 1)

        else:
            nst = int(sections['Storey'].max())
            nbays = int(sections['Bay'].max() - 1)

        # Heights and widths of the structure
        heights = np.array([0])
        widths = np.array([0])

        # Disaggregate beam and column sections
        beams = sections[(sections['Element'] == 'Beam')]
        columns = sections[(sections['Element'] == 'Column')]

        # Bottom and top nodes for the recorders
        bnode = []
        tnode = []
        for st in range(nst):
            heights = np.append(heights, heights[st] +
                                columns[(columns['Storey'] == st + 1)].iloc[0]['length'])

            if self.hingeModel == 'haselton':
                if st == 0:
                    bnode.append(int(f"{st + 1}{nbays + 1}00"))
                else:
                    bnode.append(int(f"{st + 1}{nbays + 1}20"))
                tnode.append(int(f"{st + 2}{nbays + 1}20"))

            elif self.hingeModel == 'hysteretic':
                if self.flag3d:
                    if direction == 0:
                        bnode.append(int(f"{nbays + 1}1{st}"))
                        tnode.append(int(f"{nbays + 1}1{st + 1}"))
                    else:
                        bnode.append(int(f"1{nbays + 1}{st}"))
                        tnode.append(int(f"1{nbays + 1}{st + 1}"))
                else:
                    bnode.append(int(f"{nbays + 1}{st}"))
                    tnode.append(int(f"{nbays + 1}{st + 1}"))

            else:
                raise ValueError('[EXCEPTION] Wrong lumped hinge model (should be Hysteretic or Haselton)')

        for bay in range(nbays):
            widths = np.append(widths, widths[bay] +
                               beams[(beams['Bay'] == bay + 1)].iloc[0]['length'])

        return nst, nbays, heights, widths, beams, columns, bnode, tnode

    def define_nodes(self):
        """
        defines nodes
        :return: DataFrame                      Node IDs and coordinates
        """
        # y coordinate used only for 3D modelling
        df = {'Node id': [], 'x': [], 'y': [], 'z': []}

        if self.flag3d:
            for st in range(self.nst + 1):
                for x in range(self.nbays[0] + 1):
                    for y in range(self.nbays[1] + 1):
                        df["Node id"].append(f"{x + 1}{y + 1}{st}")
                        df["z"].append(self.heights[st])
                        df["x"].append(self.widths[0][x])
                        df["y"].append(self.widths[1][y])

        else:
            for st in range(self.nst + 1):
                for bay in range(self.nbays + 1):
                    if self.hingeModel == 'haselton':
                        # 4 nodes per panel zone
                        if st == 0:
                            df['Node id'].append(f"{st + 1}{bay + 1}00")
                            df['z'].append(self.heights[st])
                            df['x'].append(self.widths[bay])
                            df['Node id'].append(f"{st + 1}{bay + 1}000")
                            df['z'].append(self.heights[st])
                            df['x'].append(self.widths[bay])

                        else:
                            beam_heights = np.array(self.beams['h'])
                            col_ext_heights = np.array(self.columns['h'][(self.columns['Position'] == 'External')])
                            if self.nbays > 1:
                                try:
                                    col_int_heights = np.array(self.columns['h'][(self.columns['Position'] == 'Internal')])
                                    df['Node id'].append(f"{st + 1}{bay + 1}10")
                                    df['z'].append(self.heights[st] - beam_heights[(st - 1)] / 2)
                                    df['x'].append(self.widths[bay])
                                    df['Node id'].append(f"{st + 1}{bay + 1}20")
                                    df['z'].append(self.heights[st])
                                    if bay == 0 or bay == self.nbays:
                                        df['x'].append(self.widths[bay] + col_ext_heights[(st - 1)] / 2)
                                    else:
                                        df['x'].append(self.widths[bay] + col_int_heights[(st - 1)] / 2)

                                    df['Node id'].append(f"{st + 1}{bay + 1}30")
                                    df['z'].append(self.heights[st] + beam_heights[(st - 1)] / 2)
                                    df['x'].append(self.widths[bay])
                                    df['Node id'].append(f"{st + 1}{bay + 1}40")
                                    df['z'].append(self.heights[st])
                                    if bay == 0 or bay == self.nbays:
                                        df['x'].append(self.widths[bay] - col_ext_heights[(st - 1)] / 2)
                                    else:
                                        df['x'].append(self.widths[bay] - col_int_heights[(st - 1)] / 2)

                                except:
                                    raise ValueError('[EXCEPTION] Internal column cross-sections not provided!')

                    else:
                        df['Node id'].append(f"{bay + 1}{st}")
                        df['z'].append(self.heights[st])
                        df['x'].append(self.widths[bay])
                        # As only a 2D frame is being considered
                        df['y'].append(0.0)

        df = pd.DataFrame.from_dict(df)
        return df
