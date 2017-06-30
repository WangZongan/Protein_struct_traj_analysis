##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2013 Stanford University and the Authors
#
# Authors: Robert McGibbon
# Contributors: Wang Zongan
#
# This program is adapted from the hbond.py of MDTraj.
# The purpose of the modification and addition is to allow user-supplied
# distance and angle cutoffs in the hbond calculation function.
# Moreover, I want to implement an additional hbond function used by VMD.
##############################################################################

##############################################################################
# Imports
##############################################################################

from __future__ import print_function, division
import numpy as np
from mdtraj.utils    import ensure_type
from mdtraj.geometry import compute_distances, compute_angles
from mdtraj.geometry import _geometry

__all__ = ['baker_hubbard_wz', 'vmd_hbond_wz']

##############################################################################
# Functions
##############################################################################

def baker_hubbard(traj,
                  distance_cutoff=0.25, angle_cutoff=2.094395,
                  freq=0.1, exclude_water=True, periodic=True, sidechain_only=False):
    """
    Identify hydrogen bonds based on cutoffs for the Donor-H...Acceptor distance and angle.

    The criterion employed is angle(DHA) > 120 and distance(HA) < 2.5 A.

    When donor the donor is 'N' and the acceptor is 'O', this corresponds to
    the definition established in [1]_. The donors considered by this method
    are NH and OH, and the acceptors considered are O and N.
    Parameters
    ----------
    traj : md.Trajectory
        An mdtraj trajectory. It must contain topology information.

    distance_cutoff: float, default=0.25 (nanometers),
        Cutoff for the distance between H and the acceptor.
    angle_cutoff: float, default=2.094395 (2.0*np.pi/3.0, radians),
        Cutoff for the angle of DHA.

    freq : float, default=0.1
        Return only hydrogen bonds that occur in greater this fraction of the
        frames in the trajectory.
    exclude_water : bool, default=True
        Exclude solvent molecules from consideration
    periodic : bool, default=True
        Set to True to calculate displacements and angles across periodic box boundaries.
    sidechain_only : bool, default=False
        Set to True to only consider sidechain-sidechain interactions.

    Returns
    -------
    hbonds : np.array, shape=[n_hbonds, 3], dtype=int
        An array containing the indices atoms involved in each of the identified
        hydrogen bonds. Each row contains three integer indices, `(d_i, h_i,
        a_i)`, such that `d_i` is the index of the donor atom, `h_i` the index
        of the hydrogen atom, and `a_i` the index of the acceptor atom involved
        in a hydrogen bond which occurs (according to the definition above) in
        proportion greater than `freq` of the trajectory.

    Notes
    -----
    Each hydrogen bond is distinguished for the purpose of this function by the
    indices of the donor, hydrogen, and acceptor atoms. This means that, for
    example, when an ARG sidechain makes a hydrogen bond with its NH2 group,
    you might see what appear like double counting of the h-bonds, since the
    hydrogen bond formed via the H_1 and H_2 are counted separately, despite
    their "chemical indistinguishably"
    Examples
    --------
    >>> md.baker_hubbard_wz(t)
    array([[  0,  10,   8],
           [  0,  11,   7],
           [ 69,  73,  54],
           [ 76,  82,  65],
           [119, 131,  89],
           [140, 148, 265],
           [166, 177, 122],
           [181, 188, 231]])
    >>> label = lambda hbond : '%s -- %s' % (t.topology.atom(hbond[0]), t.topology.atom(hbond[2]))
    >>> for hbond in hbonds:
    >>>     print label(hbond)
    GLU1-N -- GLU1-OE2
    GLU1-N -- GLU1-OE1
    GLY6-N -- SER4-O
    CYS7-N -- GLY5-O
    TYR11-N -- VAL8-O
    MET12-N -- LYS20-O
    See Also
    --------
    kabsch_sander
    References
    ----------
    .. [1] Baker, E. N., and R. E. Hubbard. "Hydrogen bonding in globular
        proteins." Progress in Biophysics and Molecular Biology
        44.2 (1984): 97-179.
    """

    if traj.topology is None:
        raise ValueError('baker_hubbard_wz requires that traj contain topology information')

    # Get the possible donor-hydrogen...acceptor triplets
    bond_triplets = _get_bond_triplets(traj.topology,
        exclude_water=exclude_water, sidechain_only=sidechain_only)

    mask, distances, angles = _compute_bounded_geometry(traj, bond_triplets,
        distance_cutoff, [1, 2], [0, 1, 2], freq=freq, periodic=periodic)

    # Find triplets that meet the criteria
    presence = np.logical_and(distances < distance_cutoff, angles > angle_cutoff)
    mask[mask] = np.mean(presence, axis=0) > freq

    return bond_triplets.compress(mask, axis=0)


def vmd_hbond(traj,
              distance_cutoff=0.35, angle_cutoff=2.094395,
              freq=0.1, exclude_water=True, periodic=True, sidechain_only=False):
    """
    Identify hydrogen bonds based on cutoffs for the Donor-H...Acceptor distance and angle.

    The criterion employed is angle(DHA) > 120 and distance(DA) < 3.5 A.

    When donor the donor is 'N' and the acceptor is 'O', this corresponds to
    the definition established in [1]_. The donors considered by this method
    are NH and OH, and the acceptors considered are O and N.
    Parameters
    ----------
    traj : md.Trajectory
        An mdtraj trajectory. It must contain topology information.

    distance_cutoff: float, default=0.35 (nanometers),
        Cutoff for the distance between donor and the acceptor.
    angle_cutoff: float, default=2.094395 (2.0*np.pi/3.0, radians),
        Cutoff for the angle of DHA.

    freq : float, default=0.1
        Return only hydrogen bonds that occur in greater this fraction of the
        frames in the trajectory.
    exclude_water : bool, default=True
        Exclude solvent molecules from consideration
    periodic : bool, default=True
        Set to True to calculate displacements and angles across periodic box boundaries.
    sidechain_only : bool, default=False
        Set to True to only consider sidechain-sidechain interactions.

    Returns
    -------
    hbonds : np.array, shape=[n_hbonds, 3], dtype=int
        An array containing the indices atoms involved in each of the identified
        hydrogen bonds. Each row contains three integer indices, `(d_i, h_i,
        a_i)`, such that `d_i` is the index of the donor atom, `h_i` the index
        of the hydrogen atom, and `a_i` the index of the acceptor atom involved
        in a hydrogen bond which occurs (according to the definition above) in
        proportion greater than `freq` of the trajectory.

    Notes
    -----
    Each hydrogen bond is distinguished for the purpose of this function by the
    indices of the donor, hydrogen, and acceptor atoms. This means that, for
    example, when an ARG sidechain makes a hydrogen bond with its NH2 group,
    you might see what appear like double counting of the h-bonds, since the
    hydrogen bond formed via the H_1 and H_2 are counted separately, despite
    their "chemical indistinguishably"
    Examples
    --------
    >>> md.vmd_hbond_wz(t)
    array([[  0,  10,   8],
           [  0,  11,   7],
           [ 69,  73,  54],
           [ 76,  82,  65],
           [119, 131,  89],
           [140, 148, 265],
           [166, 177, 122],
           [181, 188, 231]])
    """

    if traj.topology is None:
        raise ValueError('vmd_hbond_wz requires that traj contain topology information')

    # Get the possible donor-hydrogen...acceptor triplets
    bond_triplets = _get_bond_triplets(traj.topology,
        exclude_water=exclude_water, sidechain_only=sidechain_only)

    mask, distances, angles = _compute_bounded_geometry(traj, bond_triplets,
        distance_cutoff, [0, 2], [0, 1, 2], freq=freq, periodic=periodic)

    # Find triplets that meet the criteria
    presence = np.logical_and(distances < distance_cutoff, angles > angle_cutoff)
    mask[mask] = np.mean(presence, axis=0) > freq

    return bond_triplets.compress(mask, axis=0)


def _get_bond_triplets(topology, exclude_water=True, sidechain_only=False):
    '''
    Return: bond_triplets, shape (ntriplets, 3)
    Each of the elements in bond_triplets is the set of indices of (donor, H, acceptor).
    '''

    def can_participate(atom):
        # Filter waters
        if exclude_water and atom.residue.is_water:
            return False
        # Filter non-sidechain atoms
        if sidechain_only and not atom.is_sidechain:
            return False
        # Otherwise, accept it
        return True

    def get_donors(e0, e1):
        '''
        1. If get_donors('N', 'H'),

        'atoms' looks like:
        [(ASP20-N, ASP20-H), (ASP20-N, ASP20-H2), (ASP20-N, ASP20-H3), (VAL21-N, VAL21-H), ... ];

        'indices' looks like:
        [(0, 1), (0, 2), (0, 3), (14, 15), ... ].

        2. If get_donors('O', 'H'),

        'atoms' looks like:
        [(SER32-OG, SER32-HG), (TYR33-OH, TYR33-HH), (SER38-OG, SER38-HG), ... ]

        'indices' looks like:
        [(192, 193), (209, 210), (282, 283), ... ]
        '''

        # Find all matching bonds
        elems = set((e0, e1))
        atoms = [(one, two) for one, two in topology.bonds
            if set((one.element.symbol, two.element.symbol)) == elems]

        # Filter non-participating atoms
        atoms = [atom for atom in atoms
            if can_participate(atom[0]) and can_participate(atom[1])]

        # Get indices for the remaining atoms
        indices = []
        for a0, a1 in atoms:
            pair = (a0.index, a1.index)
            # make sure to get the pair in the right order, so that the index
            # for e0 comes before e1
            if a0.element.symbol == e1:
                pair = pair[::-1]
            indices.append(pair)

        return indices

    nh_donors = get_donors('N', 'H')
    oh_donors = get_donors('O', 'H')
    xh_donors = np.array(nh_donors + oh_donors) # simply concatenation of 2 lists

    if len(xh_donors) == 0:
        # if there are no hydrogens or protein in the trajectory, we get
        # no possible pairs and return nothing
        return np.zeros((0, 3), dtype=int)

    acceptor_elements = frozenset(('O', 'N'))
    acceptors = [a.index for a in topology.atoms
        if a.element.symbol in acceptor_elements and can_participate(a)]

    # Make acceptors a 2-D numpy array, shape (nacceptors, ) --> (nacceptors, 1)
    acceptors = np.array(acceptors)[:, np.newaxis]

    '''
    # Generate the cartesian product of the donors and acceptors

    >>> x = np.array([[1,2],[3,4]])
    >>> np.repeat(x, 3, axis=1)
    array([[1, 1, 1, 2, 2, 2],
           [3, 3, 3, 4, 4, 4]])
    >>> np.repeat(x, 3, axis=0)
    array([[1, 2],
           [1, 2],
           [1, 2],
           [3, 4],
           [3, 4],
           [3, 4]])

    >>> b = np.array([[1, 2], [3, 4]])
    >>> np.tile(b, (2, 1))
    array([[1, 2],
           [3, 4],
           [1, 2],
           [3, 4]])
    '''
    xh_donors_repeated = np.repeat(xh_donors, acceptors.shape[0], axis=0)
    acceptors_tiled    = np.tile(acceptors, (xh_donors.shape[0], 1))
    bond_triplets      = np.hstack((xh_donors_repeated, acceptors_tiled))

    # Filter out self-bonds
    self_bond_mask = (bond_triplets[:, 0] == bond_triplets[:, 2])
    return bond_triplets[np.logical_not(self_bond_mask), :]


def _compute_bounded_geometry(traj, triplets, distance_cutoff, distance_indices,
                              angle_indices, freq=0.0, periodic=True):
    """
    Returns a tuple include
    (1) the mask for triplets that fulfill the distance criteria frequently enough,
    (2) the actual distances calculated, and
    (3) the angles between the triplets specified by angle_indices.
    """
    # First we calculate the requested distances
    distances = compute_distances(traj, triplets[:, distance_indices], periodic=periodic)

    # Now we discover which triplets meet the distance cutoff often enough
    prevalence = np.mean(distances < distance_cutoff, axis=0)
    mask = prevalence > freq

    # Update data structures to ignore anything that isn't possible anymore
    triplets = triplets.compress(mask, axis=0)
    distances = distances.compress(mask, axis=1)

    '''
    # Calculate angles using the law of cosines

    If angle_indices = [0,1,2],
    abc_pairs = [(0, 1), (1, 2), (2, 0)].
    '''

    abc_pairs = zip(angle_indices, angle_indices[1:] + angle_indices[:1])
    abc_distances = []

    # Calculate distances (if necessary)
    for abc_pair in abc_pairs:
        if set(abc_pair) == set(distance_indices):
            abc_distances.append(distances)
        else:
            abc_distances.append(compute_distances(traj, triplets[:, abc_pair], periodic=periodic))

    # Law of cosines calculation
    a, b, c = abc_distances
    cosines = (a ** 2 + b ** 2 - c ** 2) / (2 * a * b)
    np.clip(cosines, -1, 1, out=cosines) # avoid NaN error
    angles = np.arccos(cosines)

    return mask, distances, angles
