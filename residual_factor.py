import reciprocalspaceship as rs
import numpy as np
import sys, os, getopt

# parse cmdline args in a C-like manner
try:
    opts, args = getopt.getopt(sys.argv[1:], "vDd:w:r:")
except getopt.GetoptError as err:
    print(err)
    sys.exit(1)

verbose = False
dryrun = False
dirname = ""
accumulating_avg = False
accumulating_window = 10 # frames
ref_name = ""

for o, a in opts:
    if o == "-v":
        verbose = True
    elif o == "-D":
        dryrun = True
    elif o == "-d":
        dirname = a
    elif o == "-w":
        accumulating_avg = True
        accumulating_window = int(a)
    elif o == "-r":
        ref_name = a

if dirname != "":
    try:
        os.mkdir(dirname)
    except:
        pass # ignore err due to existing dir
    os.chdir(dirname)

# load reference structural factors
ref_dataset = rs.read_mtz(f"{ref_name}.mtz")
