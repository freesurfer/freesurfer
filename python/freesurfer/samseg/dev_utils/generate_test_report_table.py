'''
This file generates a results table from a test_results json.
'''

import json
from terminaltables import AsciiTable
table_data = []
results_json_path = "/Users/ys/work/freesurfer/test_results.json"
results_dict = json.load(open(results_json_path))
for case_number, result_dict in results_dict.items():
    if 'gold versus matlab' in result_dict:
        matlab_dice = result_dict['gold versus matlab']['Dice']
        matlab_str = '{:.4f}'.format(matlab_dice)
        if matlab_dice < 0.99:
            matlab_str = '{:.4f} *'.format(matlab_dice)
    else:
        matlab_str = ""
    python_dice = result_dict['gold versus python']['Dice']
    if python_dice < 0.99:
        python_str = '{:.4f} *'.format(python_dice)
    else:
        python_str = '{:.4f}'.format(python_dice)

    table_data.append([case_number,
                       python_str,
                       matlab_str])

table_data = sorted(table_data)
header = ["ID", "Dice\nPython", "Dice\nMATLAB"]
table_data.insert(0, header)
print(AsciiTable(table_data).table)
