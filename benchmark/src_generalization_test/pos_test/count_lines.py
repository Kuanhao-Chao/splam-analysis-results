import linecache
import itertools

types = ['noN', 'N']
versions = [1,2,3,4,5]
dbs = ['GRCm39', 'Mmul_10', 'NHGRI_mPanTro3', 'TAIR10']
file_paths = [f'./5_output/spliceai_result_{version}/{db}/spliceai_all_seq.name.{type}.{db}.tsv' for db, version, type in itertools.product(dbs, versions, types)]

for file_path in file_paths:
    line_count = sum(1 for line in linecache.getlines(file_path))
    print(f"Number of lines in '{file_path}': {line_count}")