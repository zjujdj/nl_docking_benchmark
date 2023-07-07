# /home/dejun/workspace/NLDock/pymol/bin/python

# /home/dejun/workspace/NLDock/RLDOCK-master/RLDOCK -i /home/dejun/workspace/NLDock/docking/protein/3DJ0_protein.mol2 -l /home/dejun/workspace/NLDock/docking/ligand/mol2/3DJ0_ligand.mol2 -o /home/dejun/workspace/NLDock/result/docking_ran/rldock/3DJ0/3DJ0 -c 100 -n 1 -s /home/dejun/workspace/NLDock/RLDOCK-master/src/sphere.dat -r /home/dejun/workspace/NLDock/docking/ligand/mol2/3DJ0_ligand.mol2
# /home/dejun/workspace/NLDock/RLDOCK-master/RLDOCK -i /home/dejun/workspace/NLDock/docking/protein/6CC3_protein.mol2 -l ./6CC3_ligand-out.mol2 -o ./6CC3 -c 100 -n 5 -s /home/dejun/workspace/NLDock/RLDOCK-master/src/sphere.dat -r /home/dejun/workspace/NLDock/docking/ligand/mol2/6CC3_ligand.mol2
import pandas as pd
import numpy as np
import os
import multiprocessing as mp

rldock_path = '/home/dejun/workspace/NLDock/RLDOCK-master'
s_dat_path = ''
num_pose = 100


def get_docking_csv(mol2_file):
    energies, rmsds, names = [], [], []
    with open(mol2_file, 'r') as f:
        lines = f.readlines()

    for idx, line in enumerate(lines):
        if line.startswith('# Rank:'):
            names.append('Rank%s' % line.split()[2])
            rmsds.append(float(line.split()[-1]))
            energies.append(float(lines[idx + 2].split()[-1]))
    data_pd = pd.DataFrame({'name': names, 'rmsd': rmsds, 'energy': energies})
    data_pd.to_csv(mol2_file.replace('cluster.mol2', 'docking.csv'), index=False)


def rldocking(pdbid):
    try:
        print('start*************************%s***************************' % pdbid)
        protein_file = '/home/dejun/workspace/NLDock/docking/protein/%s_protein.mol2' % pdbid
        ref_file = '/home/dejun/workspace/NLDock/docking/ligand/mol2/%s_ligand.mol2' % pdbid

        dock_file = '/home/dejun/workspace/NLDock/docking/ligand_random/mol2/%s_ligand_random.mol2' % pdbid  # 对接文件
        o_prefix = '/home/dejun/workspace/NLDock/result/docking_ran/rldock/%s/%s' % (pdbid, pdbid)
        clu_mol2 = '/home/dejun/workspace/NLDock/result/docking_ran/rldock/%s/%s_cluster.mol2' % (pdbid, pdbid)

        # 清空输出路径
        res_path = '/home/dejun/workspace/NLDock/result/docking_ran/rldock/%s' % pdbid
        os.system('rm -rf %s' % res_path)
        # 创建新目录
        os.system('mkdir -p %s' % res_path)

        cmdline = '%s/RLDOCK -i %s -l %s -o %s -c %s -n 1 -s %s/src/sphere.dat -r %s' % (
        rldock_path, protein_file, dock_file, o_prefix, num_pose, rldock_path, ref_file)
        os.system(cmdline)

        if os.path.exists(clu_mol2):
            get_docking_csv(clu_mol2)
        print('start*************************%s***************************' % pdbid)
    except:
        print('*************************big error: %s***************************' % pdbid)


if __name__ == '__main__':
    import time
    st = time.time()

    protein_names = os.listdir('/home/dejun/workspace/NLDock/docking/protein/')
    protein_names = [name for name in protein_names if name.endswith('_protein.pdb')]
    pdbids = [protein_names[i].split('_')[0] for i in range(len(protein_names))]

    pdbids_ = []
    for pdbid in pdbids:
        if os.path.exists('/home/dejun/workspace/NLDock/result/docking_ran/rldock/%s/%s_docking.csv' % (pdbid, pdbid)):
            pass
        else:  # 没有成功对接的蛋白
            # 先清空对接目录
            cmdline = 'rm -rf /home/dejun/workspace/NLDock/result/docking_ran/rldock/%s' % pdbid
            os.system(cmdline)
            pdbids_.append(pdbid)

    pool_num = mp.cpu_count()
    pool = mp.Pool(10)
    pool.starmap_async(rldocking, zip(pdbids[:]))
    pool.close()
    pool.join()

    # rldocking('7TZT')

    end = time.time()
    print('total elapsed time:', end - st, 'S')