from typing import List
import gemmi
import mdtraj
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv, norm, det
DEN = 24

def real_space_transform(xyz: np.ndarray, op: gemmi.Op, unitcell_vecs: np.ndarray, wrap: bool = True, fractional: bool = False) -> np.ndarray:
    """
    Transforms an array of real coordinates by a symmetry operator and returns the new array.
    :param xyz: n * 3 array that stores the real space coordinates to be transformed
    :param op: the symop to be applied
    :param unitcell_vecs: 3 * 3 array with rows as unit cell vectors
    :param wrap: if True, fractional coordinates will be wrapped to be inside the unit box
    :param fractional: if True, returns fractional coordinate instead
    :return: n * 3 array containing real coordinates as transformed by op.
    """
    tran = np.expand_dims(np.array(op.tran) / DEN, axis=0)
    rot = np.expand_dims(inv(np.expand_dims(np.array(op.rot), axis=0)[0] / DEN), axis=0)
    deorthmat = np.array(unitcell_vecs)
    orthmat = np.expand_dims(inv(deorthmat), axis=0)
    if fractional is True:
        deorthmat = np.expand_dims(np.identity(3), axis=0)
    else:
        deorthmat = np.expand_dims(deorthmat, axis=0)
    if wrap is True:
        return ((xyz @ orthmat) @ rot + tran) % 1 @ deorthmat
    else:
        return ((xyz @ orthmat) @ rot + tran) @ deorthmat

def trajectory_transform(traj: mdtraj.Trajectory, op: gemmi.Op) -> None:
    """
    Performs an in-place transformation applied to all atomic coordinates in each frame of the input trajectory according to symmetry operator `op`
    Assumes constant unit cell vector for speed.
    Assumes CoM does't jump around near the boundary of the cell (use smart wrap on the trajectory first if needed)
    :param traj: Input trajectory that the transformation will be applied to, frame by frame.
    :param op: Symmetry operator to be applied, contains translation and rotation parts.
    :return: None.

    """
    tran = np.expand_dims(np.array(op.tran) / DEN, axis = 0)
    rot = np.expand_dims(inv(np.expand_dims(np.array(op.rot), axis=0)[0] / DEN), axis = 0)
    deorthmat = traj.unitcell_vectors[0]
    orthmat = np.expand_dims(inv(deorthmat), axis = 0)
    deorthmat = np.expand_dims(deorthmat, axis = 0)
    com = mdtraj.compute_center_of_mass(traj[0])[0:1]
    frac_com = ((com @ orthmat) @ rot + tran)
    com_shift = -(frac_com // 1)
    traj.xyz = ((traj.xyz @ orthmat) @ rot + tran + com_shift) @ deorthmat

def trajectory_revert_to_asu(traj: mdtraj.Trajectory, sg: int, chain_ids: List[int], threshold = 0.05):
    """
    Performs an in-place transformation applied to all atomic coordinates in each frame of the input trajectory to restore each chain back to the first chain (ASU)
    Assumes that each chain doesn't move too much in the unit cell so the same symop can be applied across all frames.
    Assumes reference chain is the one with id = 0.
    :param traj: Input trajectory that the transformation will be applied to.
    :return: None
    """
    frame0 =  traj[0]
    top = traj.topology
    ucv = traj.unitcell_vectors[0]
    char_length = det(ucv)**(1/3)
    ops = list(gemmi.find_spacegroup_by_number(sg).operations().sym_ops)
    ops_inv = [op.inverse() for op in ops]
    ref_sel = top.select(f"chainid 0")
    ref_com = mdtraj.compute_center_of_mass(frame0.atom_slice(ref_sel))[0]
    new_traj = traj.atom_slice(ref_sel)
    for id in chain_ids:
        if id == 0: # reference
            continue
        top_sel = top.select(f"chainid {id}")
        com = mdtraj.compute_center_of_mass(frame0.atom_slice(top_sel))[0] # (3,)
        subtraj = traj.atom_slice(top_sel)
        #print(f"Looking at chain {id}, com {com}, refcom {ref_com}")
        for op, op_inv in zip(ops, ops_inv):
            new_com = real_space_transform(com, op_inv, ucv)
            ratio =  norm(new_com - ref_com) / char_length
            #print(f"Looking at symop {op.triplet()}, ratio {ratio}")
            if ratio < threshold:
                trajectory_transform(subtraj, op_inv)
                print(f"Chain {id} corresponds to symop {op.triplet()}")
                break
        new_traj = new_traj.stack(subtraj, keep_resSeq = False)
    return new_traj

def get_closest_factors(input: int) -> (int, int):
    test_num = int(np.sqrt(input))
    while input % test_num != 0:
        test_num -= 1
    return (test_num, input // test_num)


# Courtesy of Jack's code
def compare_symops(mtz, n, original_sg=213):
    """
    The function takes a number and returns a dataset that has the x,y,z
    reflection (reciprocal space ASU) joined with the corresponding
    symmetry-related one under the given symop.

    :param mtz: A reciprocalspaceship Dataset
    :param n: the index of symop belonging to a particular space group, refer to ITC for their definitions
    :param original_sg: The space group that the dataset is assumed to come from
    :return: a Dataset including all pairs of reflections related by the specified symop
    """
    ds = mtz.copy()

    # Map P1 reflections to original ASU
    ds.spacegroup = original_sg
    asu = ds.hkl_to_asu()

    # Find common reflections
    common = asu.loc[asu["M/ISYM"]  == 1].index.intersection(asu.loc[asu["M/ISYM"]  == n].index).sort_values()
    asu = asu.loc[common]
    asu1 = asu.loc[asu["M/ISYM"] == 1].reset_index()
    asu2 = asu.loc[asu["M/ISYM"] == n].reset_index()

    return asu1.merge(asu2, on=["H", "K", "L"], suffixes=(1, 2))

# Jack Greisman's script from https://github.com/Hekstra-Lab/mdtools/blob/1f71b90d8a80d6a9d216d6a980c05221998f428a/mdtools/analysis/latticemdtrajectory.py#L37
def smartWrapMolecule(traj: mdtraj.Trajectory, indices: list) -> None:
    """
    This function applies periodic wrapping to a given set of atomic
    indices to prevent their center of mass from jumping by a unit
    cell length. Currently, it is assumed that the indices
    correspond to a molecule -- meaning a set of atoms connected by
    bonds.
    Parameters
    ----------
    indices : list of ints
        Atomic indices of positions that should be wrapped together
    """

    # Compute geometric center of coordinates
    coms = traj.xyz[:, indices, :].mean(axis=1)

    # Compute mask for integer unitcell adjustments
    mask = np.zeros(shape=(traj.n_frames, 3))

    # X-axis
    x = traj.unitcell_lengths[0, 0]
    mask[np.where(coms[:, 0] - coms[0, 0] < -1 * x / 2)[0], 0] = 1
    mask[np.where(coms[:, 0] - coms[0, 0] > x / 2)[0], 0] = -1

    # Y-axis
    y = traj.unitcell_lengths[0, 1]
    mask[np.where(coms[:, 1] - coms[0, 1] < -1 * y / 2)[0], 1] = 1
    mask[np.where(coms[:, 1] - coms[0, 1] > y / 2)[0], 1] = -1

    # Z-axis
    z = traj.unitcell_lengths[0, 2]
    mask[np.where(coms[:, 2] - coms[0, 2] < -1 * z / 2)[0], 2] = 1
    mask[np.where(coms[:, 2] - coms[0, 2] > z / 2)[0], 2] = -1

    # Update trajectory coordinates
    traj.xyz[:, indices, :] += (mask * traj.unitcell_lengths).reshape(-1, 1, 3)


def smartWrapProtein(traj: mdtraj.Trajectory, ref: mdtraj.Trajectory) -> None:
    """
    Apply smart wrapping independently to each protein molecule in
    the MD system. For now, this method identifies proteins as
    molecules with more than 100 atoms
    """
    traj = ref + traj
    for mol in traj.topology.find_molecules():
        if len(mol) > 100:
            indices = [atom.index for atom in mol]
            smartWrapMolecule(traj, indices)
    return traj[1:]


# 24 chains from the P 41 3 2 space group, we will calculate RMSD for each
def compute_RMSD_by_chain(target_traj, top, rule, chain_id_list=None):
    chains_rmsd = []
    n_frames = target_traj.n_frames
    if chain_id_list == None:
        chain_id_list = range(top.n_chains)
    for chain_id in chain_id_list:
        atoms_selection = top.select("chainid " + str(chain_id) + " and " + rule)
        if (len(atoms_selection) == 0):  # no valid atoms bound to a.a. residues
            continue
        chain_rmsd = mdtraj.rmsd(target_traj, target_traj, 0, atoms_selection)
        chains_rmsd.append(chain_rmsd)

    chains_rmsd = np.array(chains_rmsd)
    return chains_rmsd, np.mean(chains_rmsd, axis=0), np.std(chains_rmsd, axis=0)


def plot_RMSD_by_chain(chains_rmsd, rule, n_frames):
    plt.figure(figsize=(18, 12))
    plt.rc('font', size=16)
    plt.title("PDZ domain RMSD by chain (rule: {})".format(rule))
    plt.xlabel("Time (ns)")
    plt.xticks(np.arange(0, n_frames, 100), np.arange(0, n_frames,100) / 10)
    plt.ylabel("RMSD (nm)")
    for chid, chain_rmsd in enumerate(chains_rmsd):
        plt.plot(chain_rmsd, label = "chain " + str(chid))

    plt.legend(loc = 'upper right')


def plot_RMSD_by_chain_stat(avg_rmsd, std_rmsd, rule, n_frames):
    plt.figure(figsize=(18, 12))
    plt.rc('font', size=16)
    plt.title("PDZ domain RMSD by chain (rule: {})".format(rule))
    plt.xlabel("Time (ns)")
    plt.xticks(np.arange(0, n_frames, 100), np.arange(0, n_frames,100) / 10)
    plt.ylabel("RMSD (nm)")
    plt.plot(avg_rmsd, label = "averaged")
    plt.plot(avg_rmsd + std_rmsd, label = "+1 sd")
    plt.plot(avg_rmsd - std_rmsd, label = "-1 sd")
    plt.legend(loc = 'upper right')
