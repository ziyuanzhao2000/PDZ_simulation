import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys, os, getopt

def get_closest_factors(input: int) -> (int, int):
    test_num = int(np.sqrt(input))
    while input % test_num != 0:
        test_num -= 1
    return (test_num, input / test_num)

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
    df_var_diff_in_mag = pd.read_pickle(f"{input_name}_var_diff_in_mag.pkl")
    df_var_mag_of_cdiff = pd.read_pickle(f"{input_name}_var_mag_of_cdiff.pkl")
else:
    df_var_diff_in_mag = pd.read_pickle("var_diff_in_mag.pkl")
    df_var_mag_of_cdiff = pd.read_pickle("var_mag_of_cdiff.pkl")

dict = {"var_diff_in_mag": df_var_diff_in_mag, "var_mag_of_cdiff": df_var_mag_of_cdiff}
for key in dict.keys():
    df = dict[key]
    nx, ny = get_closest_factors(df.shape[1])
    fig, axes = plt.subplots(nx, ny, figsize=(nx * subplot_width, ny * subplot_height), tight_layout=True)
    df.plot(ax=axes, subplots=True, rot=60)

