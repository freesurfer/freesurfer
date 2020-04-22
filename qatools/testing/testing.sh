# should be run from within the 'testing' directory, otherwise adjust pathnames

python3 ../qatools.py --subjects_dir data --output_dir output/test01 --subjects 091

python3 ../qatools.py --subjects_dir data --output_dir output/test02 --subjects 129 130

python3 ../qatools.py --subjects_dir data --output_dir output/test03

python3 ../qatools.py --subjects_dir data --output_dir output/test04 --screenshots

python3 ../qatools.py --subjects_dir data --output_dir output/test05 --fornix

python3 ../qatools.py --subjects_dir data --output_dir output/test06 --outlier

python3 ../qatools.py --subjects_dir data --output_dir output/test07 --outlier --outlier-table data/buckner_AsegExampleNorms.csv

python3 ../qatools.py --subjects_dir data --output_dir output/test08 --screenshots --fornix --outlier

