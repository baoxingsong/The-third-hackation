# -*- coding: utf-8 -*-
# @Author: wjq
# @Date:   2021-07-12 15:16:31
# @Last Modified by:   wjq
# @Last Modified time: 2022-03-14 00:45:16
"""
calculata ks

"""
import os
import glob
import subprocess
from multiprocessing import Pool


class RunBioSoftware:
    def __init__(self, num_process):
        self.num_process = num_process
        self.main()

    @staticmethod
    def check_file(file_path):
        n = os.path.isfile(file_path)
        fn_name = os.path.basename(file_path)
        assert n, f"{fn_name}, File not exist!!!"

    def run_ks(self, path, blk, cds1, cds2):
        # spec1, spec2 = path.split('_')
        spec1, spec2 = path, path
        ks_sf = f'{spec1}_{spec2}.ks.txt'
        if not os.path.isdir(path):
            os.mkdir(path)
        if not os.path.isdir('ks/'):
            os.mkdir('ks/')
        cds_f = f'{path}/{spec1}_{spec2}.cds'
        if spec1 == spec2:
            os.system(f'cp {cds1} {cds_f}')
        else:
            os.system(f'cat {cds1} {cds2} > {cds_f}')

        p = subprocess.call(f'perl calculate.Ks.pl {path} {spec1} {blk} {cds_f} {ks_sf}', shell=True)
        if p != 0:
            print(f'Warning!!! {p}')
        # os.system(f'rm -rf {path}')

    def main(self):
        print('start!!!')
        p = Pool(self.num_process)
        files = glob.glob(f'./blk/*.pair')
        for fn in files:
            if 'partially' in fn or 'less' in fn:
                continue
            name = os.path.basename(fn).split('.')[0]
            # spec1, spec2 = name.split('_')
            spec1, spec2 = name, name
            cds1, cds2 = f'./cds/{spec1}.cds', f'./cds/{spec2}.cds'
            self.check_file(cds1), self.check_file(cds2) 
            print(name, fn, cds1, cds2)
            p.apply_async(self.run_ks, args=((name, fn, cds1, cds2)))


        p.close()
        p.join()


if __name__ == '__main__':

    process = 6  # default
    RunBioSoftware(process)


