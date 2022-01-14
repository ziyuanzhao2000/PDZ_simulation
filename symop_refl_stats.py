import reciprocalspaceship as rs
from reciprocalspaceship.utils import  to_structurefactor as to_sf
import numpy as np
import pandas as pd
import sys, os, getopt
import gemmi
from utils import compare_symops


# parse cmdline args in a C-like manner
try:
    opts, args = getopt.getopt(sys.argv[1:], "v:d:i:o:g:n:D:t:")
except getopt.GetoptError as err:
    print(err)
    sys.exit(1)

input_name = ""
output_name = ""
dir_name = ""
sg = 213
n_frames = 100
d_frame = 10
t_per_frame = 1 # ns
verbose = False

for o, a in opts:
    if o == "-i":
        input_name = a
    elif o == "-d":
        dir_name = a
    elif o == "-o":
        output_name = a
    elif o == "-g":
        sg = int(a)
    elif o == "-n":
        n_frames = int(a)
    elif o == "-D":
        d_frame = int(a)
    elif o == "-t":
        t_per_frame = int(a)
    elif o == "-v":
        verbose = True

try:
    dataset = rs.read_mtz(f"{input_name}.mtz")
    n_reflections = dataset.shape[0]

except:
    print("Can't load the mtz file for the designated chain at frame 0!")
    sys.exit(1)

time_series_index = list(np.arange(0,n_frames,d_frame)*t_per_frame)
for i, symop in enumerate(gemmi.find_spacegroup_by_number(sg).operations()):
    var_diff_in_mag = []
    var_mag_of_cdiff = []
    for frame in range(0, n_frames, d_frame):
        symop_dataset = compare_symops(dataset, i + 1, sg)
        var_diff_in_mag.append(np.var(np.abs(symop_dataset['FMODEL1'] - symop_dataset['FMODEL2'])))
        var_mag_of_cdiff.append(np.var(np.abs(to_sf(symop_dataset['FMODEL1'], symop_dataset['PHIFMODEL1']) -\
                                         to_sf(symop_dataset['FMODEL2'], symop_dataset['PHIFMODEL2']))))
    df1 = pd.DataFrame(var_diff_in_mag, columns=[symop.triplet], index=time_series_index)
    df2 = pd.DataFrame(var_mag_of_cdiff, columns=[symop.triplet], index=time_series_index)
    if i == 0: # initialize
        df_var_diff_in_mag = df1
        df_var_mag_of_cdiff = df2
    else: # bind
        df_var_diff_in_mag = pd.concat([df_var_diff_in_mag, df1], axis=1)
        df_var_mag_of_cdiff = pd.concat([df_var_mag_of_cdiff, df2], axis=1)
    if verbose:
        print(f"Finished processing symop no.{i+1}")

if verbose:
    print(df_var_diff_in_mag, df_var_mag_of_cdiff)

if output_name != "":
    df_var_diff_in_mag.to_pickle(f"{output_name}_var_diff_in_mag.pkl")
    df_var_mag_of_cdiff.to_pickle(f"{output_name}_var_mag_of_cdiff.pkl")
else:
    df_var_diff_in_mag.to_pickle("var_diff_in_mag.pkl")
    df_var_mag_of_cdiff.to_pickle("var_mag_of_cdiff.pkl")
