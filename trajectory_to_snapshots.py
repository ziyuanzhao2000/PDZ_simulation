import mdtraj

frame_skip = 10 # sample trajectory every 'frame_skip' frames
n_chains = 48 # 2 * 24

# Loads the .h5 trajectory file from the post-squeeze production run using OpenMM
target_traj = mdtraj.load("/Users/zhaoziyuan/Dropbox/Fall21/Research/Doeke_lab/PDZ_simulation/liganded_squeezed_production1.h5")
top = target_traj.topology

# slice
target_traj = target_traj[::frame_skip]
ref_traj = mdtraj.load("/Users/zhaoziyuan/Dropbox/Fall21/Research/Doeke_lab/PDZ_simulation/pdz3_cript_ds1_refine_193_fixed.pdb")

for i, frame in enumerate(target_traj):
    for chain_id in range(0, 48, 2):
        chain = frame.atom_slice(top.select("chainid " + str(chain_id) + " or chainid " + str(chain_id + 1) + " and protein"))
        chain.superpose(ref_traj) # align using all atoms, water, ions, and protein included
        chain.save_pdb("liganded_pdz_traj_" + str(i) + "_" + str(chain_id / 2) + ".pdb")