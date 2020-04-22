# should be run from within the 'testing' directory, otherwise adjust pathnames

from qatoolspython import qatoolspython

qatoolspython.run_qatools(subjects_dir='data', output_dir='output/test11', subjects=['091'])

qatoolspython.run_qatools(subjects_dir='data', output_dir='output/test12', subjects=['129', '130'])

qatoolspython.run_qatools(subjects_dir='data', output_dir='output/test13')

qatoolspython.run_qatools(subjects_dir='data', output_dir='output/test14', screenshots=True)

qatoolspython.run_qatools(subjects_dir='data', output_dir='output/test15', fornix=True)

qatoolspython.run_qatools(subjects_dir='data', output_dir='output/test16', outlier=True)

qatoolspython.run_qatools(subjects_dir='data', output_dir='output/test17', outlier=True, outlier_table='data/buckner_AsegExampleNorms.csv')

qatoolspython.run_qatools(subjects_dir='data', output_dir='output/test18', shape=True, screenshots=True, fornix=True, outlier=True)
