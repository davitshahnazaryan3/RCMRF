"""
Adapted from openseespy
"""

import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import openseespy.postprocessing.internal_database_functions as idbf
import openseespy.postprocessing.internal_plotting_functions as ipltf

for line in range(0, len(sys.argv)):
    if "ipykernel_launcher.py" in sys.argv[line]:
        matplotlib.use('nbagg')
        break
    else:
        pass

# All the plotting related definitions start here.
ele_style = {'color': 'black', 'linewidth': 1, 'linestyle': '-'}  # elements
node_style = {'color': 'black', 'marker': 'o', 'facecolor': 'black', 'linewidth': 0.}
node_style_animation = {'color': 'black', 'marker': 'o', 'markersize': 2., 'linewidth': 0.}

node_text_style = {'fontsize': 8, 'fontweight': 'regular', 'color': 'blue'}
ele_text_style = {'fontsize': 8, 'fontweight': 'bold', 'color': 'darkred'}

WireEle_style = {'color': 'black', 'linewidth': 1, 'linestyle': ':'}  # elements
Eig_style = {'color': 'red', 'linewidth': 1, 'linestyle': '-'}  # elements


def plot_model(*argv, Model="none"):
    """
    Command: plot_model(<"nodes">,<"elements">,<Model="ModelName">)

    nodes	: String, Optional, takes user input to show node tags on the model
    elements: String, Optional, takes user input to show element tags on the model
    Model	: Optional input for the name of the model used in createODB() to read the modeshape data from.
                  The default is "none" and the mode shape is plotted from the active model.

    Matplotlib rendering is faster when tags are not displayed.

    """
    # Default values
    show_node_tags = 'no'
    show_element_tags = 'no'

    # Process inputs to allow for backwards compatibility
    if len(argv) > 0:
        if any(nodeArg in argv for nodeArg in ["nodes", "Nodes", "node", "Node"]):
            show_node_tags = 'yes'
        if any(eleArg in argv for eleArg in ["elements", "Elements", "element", "Element"]):
            show_element_tags = 'yes'
        if show_node_tags == "no" and show_element_tags == "no":
            raise Exception(
                "Wrong input arguments. Command should be plot_model(<'node'>,<'element'>,Model='model_name')")

    # TODO make this a function?
    # Check if their is an output database or not.
    if Model == "none":
        print("No Model_ODB specified, trying to get data from the active model.")
        try:
            nodeArray, elementArray = idbf._getNodesandElements()
        except:
            raise Exception("No Model_ODB specified. No active model found.")
    else:
        print("Reading data from the " + Model + "_ODB.")
        try:
            nodeArray, elementArray = idbf._readNodesandElements(Model)
        except:
            raise Exception("No Model_ODB found. No active model found.")

    nodetags = nodeArray[:, 0]

    def nodecoords(nodetag):
        """
        Returns an array of node coordinates: works like nodeCoord() in opensees.
        """
        i, = np.where(nodeArray[:, 0] == float(nodetag))
        return nodeArray[int(i), 1:]

    # Check if the model is 2D or 3D
    if len(nodecoords(nodetags[0])) == 2:
        print('2D model')
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        for ele in elementArray:
            eleTag = int(ele[0])
            Nodes = ele[1:]

            if len(Nodes) == 2:
                # 2D beam-column elements
                iNode = nodecoords(Nodes[0])
                jNode = nodecoords(Nodes[1])

                ipltf._plotBeam2D(iNode, jNode, ax, show_element_tags, eleTag, "solid")

            if len(Nodes) == 3:
                # 2D Planer three-node shell elements
                iNode = nodecoords(Nodes[0])
                jNode = nodecoords(Nodes[1])
                kNode = nodecoords(Nodes[2])

                ipltf._plotTri2D(iNode, jNode, kNode, ax, show_element_tags, eleTag, ele_style, fillSurface='yes')

            if len(Nodes) == 4:
                # 2D Planer four-node shell elements
                iNode = nodecoords(Nodes[0])
                jNode = nodecoords(Nodes[1])
                kNode = nodecoords(Nodes[2])
                lNode = nodecoords(Nodes[3])

                ipltf._plotQuad2D(iNode, jNode, kNode, lNode, ax, show_element_tags, eleTag, ele_style,
                                  fillSurface='yes')

        if show_node_tags == 'yes':
            for node in nodetags:
                ax.text(nodecoords(node)[0] * 1.02, nodecoords(node)[1] * 1.02, str(int(node)),
                        **node_text_style)  # label nodes

            ax.scatter(nodeArray[:, 1], nodeArray[:, 2], **node_style)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')

    else:
        print('3D model')
        fig = plt.figure(figsize=(20, 14), dpi=80)
        ax = fig.add_subplot(1, 1, 1, projection='3d')

        for ele in elementArray:
            eleTag = int(ele[0])
            Nodes = ele[1:]

            if len(Nodes) == 2:
                # 3D beam-column elements
                iNode = nodecoords(Nodes[0])
                jNode = nodecoords(Nodes[1])

                ipltf._plotBeam3D(iNode, jNode, ax, show_element_tags, eleTag, "solid")

            if len(Nodes) == 4:
                # 3D four-node Quad/shell element
                iNode = nodecoords(Nodes[0])
                jNode = nodecoords(Nodes[1])
                kNode = nodecoords(Nodes[2])
                lNode = nodecoords(Nodes[3])

                ipltf._plotQuad3D(iNode, jNode, kNode, lNode, ax, show_element_tags, eleTag, ele_style,
                                  fillSurface='yes')

            if len(Nodes) == 8:
                # 3D eight-node Brick element
                # Nodes in CCW on bottom (0-3) and top (4-7) faces resp
                iNode = nodecoords(Nodes[0])
                jNode = nodecoords(Nodes[1])
                kNode = nodecoords(Nodes[2])
                lNode = nodecoords(Nodes[3])
                iiNode = nodecoords(Nodes[4])
                jjNode = nodecoords(Nodes[5])
                kkNode = nodecoords(Nodes[6])
                llNode = nodecoords(Nodes[7])

                ipltf._plotCubeVol(iNode, jNode, kNode, lNode, iiNode, jjNode, kkNode, llNode, ax, show_element_tags,
                                   eleTag, 'solid', fillSurface='yes')

        if show_node_tags == 'yes':
            for node in nodetags:
                ax.text(nodecoords(node)[0] * 1.02, nodecoords(node)[1] * 1.02, nodecoords(node)[2] * 1.02,
                        str(int(node)), **node_text_style)  # label nodes

            ax.scatter(nodeArray[:, 1], nodeArray[:, 2], nodeArray[:, 3], **node_style)  # show nodes

    ipltf._setStandardViewport(fig, ax, nodeArray[:, 1:], len(nodecoords(nodetags[0])))

    ax.set_xlabel('X, [m]')
    ax.set_ylabel('Y, [m]')
    ax.set_zlabel('Z, [m]')

    ax.set_xlim([-5, 40])
    ax.set_ylim([-5, 25])
    ax.set_zlim([0, 15])
    plt.axis('on')
    plt.show()
    return fig, ax
