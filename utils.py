import mdtraj
import numpy as np
import matplotlib.pyplot as plt


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
