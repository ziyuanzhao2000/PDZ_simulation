## Step 0. Initializing and house keeping
import os, sys, getopt
## step 1
from msmbuilder.io import gather_metadata, save_meta, NumberedRunsParser

## step 2.
import mdtraj as md
from msmbuilder.featurizer import DihedralFeaturizer
from msmbuilder.io import load_meta, preload_tops, save_trajs, save_generic
from multiprocessing import Pool

## Step 3.
from msmbuilder.io import load_trajs, save_trajs, save_generic
from msmbuilder.decomposition import tICA

## Step 4.
from msmbuilder.io import load_trajs, save_trajs, save_generic
from msmbuilder.cluster import MiniBatchKMeans

dir_name = "trajs"
traj_name = ""
top_name = "top"
delta_t = 100 # picoseconds per step
format = "xtc"
dim = 2 # For clustering tICA space
n_clusters = 5 # For clustering tICA space
phase = 0
verbose = False

# parse cmdline args in a C-like manner
try:
    opts, args = getopt.getopt(sys.argv[1:], "vDp:d:i:t:")
except getopt.GetoptError as err:
    print(err)
    sys.exit(1)

for o, a in opts:
    if o == "-v":
        verbose = True
    elif o == "-D":
        dryrun = True
    elif o == "-i":
        traj_name = a
    elif o == "-p":
        phase = int(a)
    elif o == "-d":
        dir_name = a
    elif o == "-t":
        top_name = a

if dir_name != "":
    try:
        os.mkdir(dir_name)
    except:
        pass # ignore err due to existing dir
    os.chdir(dir_name)

## Step 0. First obtain and process trajectory data
if phase == 0:
    phase = 1

## Step 1. Construct and save the dataframe
if phase == 1:
    parser = NumberedRunsParser(
        traj_fmt=f"{traj_name}_{{run}}.{format}",
        top_fn=f"{top_name}.pdb",
        step_ps=delta_t,
    )
    meta = gather_metadata(f"{traj_name}_*.{format}", parser)
    save_meta(meta)
    if verbose:
        print("Finished gathering data and constructing the metadata dataframe")
    phase = 2

## Step 2. Featurize the trajectories
if phase == 2:
    ## Load
    meta = load_meta()
    tops = preload_tops(meta)
    dihed_feat = DihedralFeaturizer()

    ## Featurize logic
    def feat(irow):
        i, row = irow
        traj = md.load(row['traj_fn'], top=tops[row['top_fn']])
        feat_traj = dihed_feat.partial_transform(traj)
        return i, feat_traj

    ## Do it in parallel
    with Pool() as pool:
        dihed_trajs = dict(pool.imap_unordered(feat, meta.iterrows()))

    ## Save
    save_trajs(dihed_trajs, 'ftrajs', meta)
    save_generic(dihed_feat, 'featurizer.pickl')

    if verbose:
        print("Finished featurizing the trajectories")
    phase = 3

## Step 3. Do tICA
if phase == 3:
    ## Load
    tica = tICA(n_components=5, lag_time=100, kinetic_mapping=True)
    meta, ftrajs = load_trajs("ftrajs")

    ## Fit
    tica.fit(ftrajs.values())

    ## Transform
    ttrajs = {}  ## ttrajs for transformed trajs
    for k, v in ftrajs.items():
        ttrajs[k] = tica.partial_transform(v)

    ## Save
    save_trajs(ttrajs, 'ttrajs', meta)
    save_generic(tica, 'tica.pickl')
    if verbose:
        print("Finished tICA projection")
    phase = 4

## Step 4. Analysis
if phase == 4:
    # Cluster the tICA space
    ## Load
    meta, ttrajs = load_trajs('ttrajs')

    ## Fit
    kmeans = MiniBatchKMeans(n_clusters=n_clusters)
    kmeans.fit([traj[:, :dim] for traj in ttrajs.values()])

    ## Transform
    ktrajs = {}
    for k, v in ttrajs.items():
        ktrajs[k] = kmeans.partial_transform(v[:, :dim])

    ## Save
    print(kmeans.summarize()) ## Not implemented??
    save_trajs(ktrajs, 'ktrajs', meta)
    save_generic(kmeans, 'kmeans.pickl')
    if verbose:
        print("Finished post-tICA analysis")
    phase = 5