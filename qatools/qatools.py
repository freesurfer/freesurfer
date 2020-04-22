#!/usr/bin/env python3

# say hello
print("")
print("-----------------------------")
print("qatools-python")
print("-----------------------------")
print("")

# imports
from qatoolspython import qatoolspython

# parse arguments
subjects_dir, output_dir, subjects, shape, screenshots, fornix, outlier, outlier_table = qatoolspython._parse_arguments()

# run qatools
qatoolspython.run_qatools(subjects_dir, output_dir, subjects, shape, screenshots, fornix, outlier, outlier_table)

