"""
Creates geometry of the model
"""
import pandas as pd
import numpy as np


class Geometry:
    def __init__(self, sections, hingeModel):
        """
        initializes geometry creation
        :param sections: DataFrame              Section properties
        :param hingeModel: str                  Hinge model type (Haselton (4 nodes) or Hysteretic (2 nodes))
        """
        self.sections = sections
        self.hingeModel = hingeModel.lower()

        # Get the number of storeys and bays in the direction of seismic action
        if self.hingeModel == 'haselton':
            self.nst = max(self.sections['Storey'])
            self.nbays = int(len(self.sections[(self.sections['Element'] == 'Column')]) / self.nst - 1)
        else:
            self.nst = self.sections['Storey'].max()
            self.nbays = self.sections['Bay'].max() - 1

        # Heights and widths of the structure
        self.heights = np.array([0])
        self.widths = np.array([0])

        # Disaggregate beam and column sections
        self.beams = self.sections[(self.sections['Element'] == 'Beam')]
        self.columns = self.sections[(self.sections['Element'] == 'Column')]

        # Bottom and top nodes for the recorders
        self.bnode = []
        self.tnode = []
        for st in range(self.nst):
            self.heights = np.append(self.heights, self.heights[st] +
                                     self.columns[(self.columns['Storey'] == st + 1)].iloc[0]['length'])
            if self.hingeModel == 'haselton':
                if st == 0:
                    self.bnode.append(int(f"{st + 1}{self.nbays + 1}00"))
                else:
                    self.bnode.append(int(f"{st + 1}{self.nbays + 1}20"))
                self.tnode.append(int(f"{st + 2}{self.nbays + 1}20"))
            elif self.hingeModel == 'hysteretic':
                self.bnode.append(int(f"{self.nbays + 1}{st}"))
                self.tnode.append(int(f"{self.nbays + 1}{st + 1}"))
            else:
                raise ValueError('[EXCEPTION] Wrong lumped hinge model (should be Hysteretic or Haselton)')

        for bay in range(self.nbays):
            self.widths = np.append(self.widths, self.widths[bay] +
                                    self.beams[(self.beams['Bay'] == bay + 1)].iloc[0]['length'])

    def define_nodes(self):
        """
        defines nodes
        :return: DataFrame                      Node IDs and coordinates
        """
        df = {'Node id': [], 'x': [], 'z': []}
        for st in range(self.nst + 1):
            for bay in range(self.nbays + 1):
                if self.hingeModel == 'haselton':
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

        df = pd.DataFrame.from_dict(df)
        return df
