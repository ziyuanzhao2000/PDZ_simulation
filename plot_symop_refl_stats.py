import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys, os, getopt
from utils import get_closest_factors

# parse cmdline args in a C-like manner
try:
    opts, args = getopt.getopt(sys.argv[1:], "d:i:o:")
except getopt.GetoptError as err:
    print(err)
    sys.exit(1)

input_name = ""
output_name = ""
dir_name = ""
subplot_width = 6
subplot_height = 3

for o, a in opts:
    if o == "-i":
        input_name = a
    elif o == "-d":
        dir_name = a
    elif o == "-o":
        output_name = a

if dir_name != "":
    try:
        os.mkdir(dir_name)
    except:
        pass # ignore err due to existing dir
    os.chdir(dir_name)

if input_name != "":
    df_var_diff_in_mag = pd.read_pickle(f"{input_name}_var_diff_in_mag.pkl").fillna(0)
    df_var_mag_of_cdiff = pd.read_pickle(f"{input_name}_var_mag_of_cdiff.pkl").fillna(0)
else:
    df_var_diff_in_mag = pd.read_pickle("var_diff_in_mag.pkl").fillna(0)
    df_var_mag_of_cdiff = pd.read_pickle("var_mag_of_cdiff.pkl").fillna(0)

dict = {"var_diff_in_mag": ("Magnitude of Complex Difference", df_var_diff_in_mag),
        "var_mag_of_cdiff": ("Difference of Complex Magnitudes", df_var_mag_of_cdiff)}
for key in dict.keys():
    title, df = dict[key]
    nx, ny = get_closest_factors(df.shape[1])
    print(nx, ny, df.shape)
    fig, axes = plt.subplots(nx, ny, figsize=(nx * subplot_width, ny * subplot_height), tight_layout=True,
                             sharex='col', sharey='row')
    df.plot(ax=axes, subplots=True, rot=60, title=list(df.columns),
            xlabel = "Cumulative Time (ns)", ylabel = "Variance", legend = False)
    fig.suptitle(title, fontsize=16)
    if input_name != "":
        fig.savefig(f"{output_name}_{key}.pdf")
    else:
        fig.savefig(f"{key}.pdf")