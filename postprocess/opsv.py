import openseespy.opensees as ops
import numpy as np


def section_force_diagram_2d(Ew, nep=17):
    ele_tags = ops.getEleTags()

    forces = {}
    for ele_tag in ele_tags:

        # by default no element load
        eload_data = ['', 0., 0.]
        if ele_tag in Ew:
            eload_data = Ew[ele_tag]

        nd1, nd2 = ops.eleNodes(ele_tag)

        # element x, y coordinates
        ex = np.array([ops.nodeCoord(nd1)[0],
                       ops.nodeCoord(nd2)[0]])
        ey = np.array([ops.nodeCoord(nd1)[1],
                       ops.nodeCoord(nd2)[1]])

        Lxy = np.array([ex[1]-ex[0], ey[1]-ey[0]])

        pl = ops.eleResponse(ele_tag, 'localForces')

        forces[ele_tag], _ = section_force_distribution_2d(ex, ey, pl, nep, eload_data)

    return forces


def section_force_distribution_2d(ex, ey, pl, nep=2, ele_load_data=None):
    """
    Calculate section forces (N, V, M) for an elastic 2D Euler-Bernoulli beam.

    Input:
    ex, ey - x, y element coordinates in global system
    nep - number of evaluation points, by default (2) at element ends
    ele_load_list - list of transverse and longitudinal element load
      syntax: [ele_load_type, Wy, Wx]
      For now only '-beamUniform' element load type is acceptable

    Output:
    s = [N V M]; shape: (nep,3)
        section forces at nep points along local x
    xl: coordinates of local x-axis; shape: (nep,)

    Use it with dia_sf to draw N, V, M diagrams.

    TODO: add '-beamPoint' element load type
    """

    # eload_type, Wy, Wx = ele_load_data[0], ele_load_data[1], ele_load_data[2]
    if ele_load_data is None:
        ele_load_data = ['-beamUniform', 0., 0.]

    Wy, Wx = ele_load_data[1], ele_load_data[2]

    nlf = len(pl)
    if nlf == 2:  # trusses
        N_1 = pl[0]
    elif nlf == 6:  # plane frames
        # N_1, V_1, M_1 = pl[0], pl[1], pl[2]
        N_1, V_1, M_1 = pl[:3]
    else:
        print('\nWarning! Not supported. Number of nodal forces: {nlf}')

    Lxy = np.array([ex[1]-ex[0], ey[1]-ey[0]])
    L = np.sqrt(Lxy @ Lxy)

    xl = np.linspace(0., L, nep)
    one = np.ones(nep)

    N = -1.*(N_1 * one + Wx * xl)

    if nlf == 6:
        V = V_1 * one + Wy * xl
        M = -M_1 * one + V_1 * xl + 0.5 * Wy * xl**2
        s = np.column_stack((N, V, M))
    elif nlf == 2:
        s = np.column_stack((N))

    # if eload_type == '-beamUniform':
    # else:

    return s, xl


def section_force_diagram_3d(Ew, nep=17):
    ele_tags = ops.getEleTags()

    forces = {}
    for i, ele_tag in enumerate(ele_tags):
        if ele_tag not in Ew:
            # skip the columns
            continue

        # by default no element load
        eload_data = ['-beamUniform', 0., 0., 0.]
        if ele_tag in Ew:
            eload_data = Ew[ele_tag]

        nd1, nd2 = ops.eleNodes(ele_tag)

        # element x, y coordinates
        ex = np.array([ops.nodeCoord(nd1)[0],
                       ops.nodeCoord(nd2)[0]])
        ey = np.array([ops.nodeCoord(nd1)[1],
                       ops.nodeCoord(nd2)[1]])
        ez = np.array([ops.nodeCoord(nd1)[2],
                       ops.nodeCoord(nd2)[2]])

        pl = ops.eleResponse(ele_tag, 'localForces')

        s_all, xl = section_force_distribution_3d(ex, ey, ez, pl, nep,
                                                  eload_data)
        forces[ele_tag] = s_all

    return forces


def section_force_distribution_3d(ex, ey, ez, pl, nep=2,
                                  ele_load_data=None):
    # eload_type = ele_load_data[0]
    if ele_load_data is None:
        ele_load_data = ['-beamUniform', 0., 0., 0.]

    Wy, Wz, Wx = ele_load_data[1], ele_load_data[2], ele_load_data[3]

    N1, Vy1, Vz1, T1, My1, Mz1 = pl[:6]

    Lxyz = np.array([ex[1]-ex[0], ey[1]-ey[0], ez[1]-ez[0]])
    L = np.sqrt(Lxyz @ Lxyz)

    xl = np.linspace(0., L, nep)
    one = np.ones(nep)

    N = -1.*(N1*one + Wx*xl)
    Vy = Vy1*one + Wy*xl
    Vz = Vz1*one + Wz*xl
    T = -T1*one
    Mz = -Mz1*one + Vy1*xl + 0.5*Wy*xl**2
    My = My1*one + Vz1*xl + 0.5*Wz*xl**2

    s = np.column_stack((N, Vy, Vz, T, My, Mz))

    return s, xl
