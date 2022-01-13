# This simply loads one tICA model and projects a different set of trajectories, ah, and plots for convenience

from msmbuilder.io import load_trajs, save_trajs, load_generic, save_generic
import os, sys, getopt

model_dir_name = ""
traj_dir_name = ""

# parse cmdline args in a C-like manner
try:
    opts, args = getopt.getopt(sys.argv[1:], "i:I:")
except getopt.GetoptError as err:
    print(err)
    sys.exit(1)

for o, a in opts:
    if o == "-i":
        model_dir_name = a
    elif o == "-I":
        traj_dir_name = a

tica = load_generic(f'{model_dir_name}/tica.pickl')
meta, ftrajs = load_trajs(f"{traj_dir_name}/ftrajs", meta=f"{traj_dir_name}/meta.pandas.pickl")

ttrajs = {}
for k, v in ftrajs.items():
    ttrajs[k] = tica.partial_transform(v)

save_trajs(ttrajs, f'{traj_dir_name}/{model_dir_name}_ttrajs', meta)

os.system(f"python tica-plot.py -d {traj_dir_name} -t {model_dir_name}_ttrajs")