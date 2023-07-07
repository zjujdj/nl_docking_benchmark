## !/opt/anaconda3/5.3.1/bin/python3
# -*- coding: utf-8 -*-
# author huifeng
# 用 vina/gpu-1.5.3  mgltools/1.5.7进行重对接
# py27 , module load pymol
# module purge && source activate /home/dejun/miniconda2/envs/py27 && module load pymol/2.5.4 && python autodock_vina_pipeline.py
# ========================
import csv
import os
import pymol
import multiprocessing as mp
from functools import partial
import numpy as np


def get_mol2_size_center(ligand_path):  # ligand_path是mol2格式文件
    n_atoms = 0
    x_coord = []
    y_coord = []
    z_coord = []
    with open(ligand_path, 'r') as f:
        lines = f.readlines()
        for idx, line in enumerate(lines):
            if line.startswith('@<TRIPOS>ATOM'):
                st_idx = idx
            if line.startswith('@<TRIPOS>BOND'):
                end_idx = idx
        for idx, line in enumerate(lines[st_idx + 1:end_idx]):
            n_atoms = n_atoms + 1
            # x_coord.append(float(line.split()[2]))
            # y_coord.append(float(line.split()[3]))
            # z_coord.append(float(line.split()[4]))
            x_coord.append(float(line[16:26].strip()))  # x坐标
            y_coord.append(float(line[26:36].strip()))  # y坐标
            z_coord.append(float(line[36:46].strip()))  # z坐标
    x_max, x_min = np.max(x_coord), np.min(x_coord)
    y_max, y_min = np.max(y_coord), np.min(y_coord)
    z_max, z_min = np.max(z_coord), np.min(z_coord)
    x_size = (x_max - x_min) * 5
    y_size = (y_max - y_min) * 5
    z_size = (z_max - z_min) * 5
    x_center, y_center, z_center = sum(x_coord) / n_atoms, sum(y_coord) / n_atoms, sum(z_coord) / n_atoms
    return x_size, y_size, z_size, x_center, y_center, z_center


def write_mol2(path, pdb_id, ligand_name):
    cmdline = 'cd %s/result/docking_ran/vina/%s &&' % (path, pdb_id)
    cmdline += 'module purge &&'
    cmdline += 'module load openbabel &&'
    cmdline += 'obabel -ipdbqt %s_docked.pdbqt -omol2 -O %s_docked.mol2' % (
        ligand_name.split('.')[0], ligand_name.split('.')[0])
    os.system(cmdline)

    f = open('%s/result/docking_ran/vina/%s/%s_docked.mol2' % (path, pdb_id, ligand_name.split('.')[0]), 'r')
    lines = f.read()
    f.close()
    lines = lines.split('@<TRIPOS>MOLECULE\n')
    lines = lines[1:]  # 第一个元素为''
    for i, content in enumerate(lines):
        key = '%s_ligand_vina_%s.mol2' % (pdb_id, i + 1)
        out_file = '%s/result/docking_ran/vina/%s/%s' % (path, pdb_id, key)
        with open(out_file, 'w') as f:
            f.write('@<TRIPOS>MOLECULE\n' + content)

        # 利用pymol为分子加上h进行后处理
        pymol.cmd.load(
            '%s/result/docking_ran/vina/%s/%s' % (path, pdb_id, key))
        pymol.cmd.h_add("all")
        pymol.cmd.save(
            '%s/result/docking_ran/vina/%s/%s' % (path, pdb_id, key))
        pymol.cmd.delete("all")


def vina_score(path, pdb_id, num_pose):
    ###### 输入文件以及输出路径定义
    ligand_name = str(pdb_id + '_ligand.mol2')
    protein_name = pdb_id + '_protein'
    ligand_raw = path + '/docking/ligand/mol2/' + ligand_name
    protein_raw = path + '/docking/protein/' + protein_name + '.pdb'
    
    dock_name = str(pdb_id + '_ligand_random.mol2')  # 需要对接的分子结构
    dock_raw = path + '/docking/ligand_random/mol2/' + dock_name
    
    if os.path.exists((path + '/result/docking_ran/vina/%s' % pdb_id)):  
        cmdline = 'rm -rf %s' % (path + '/result/docking_ran/vina/%s' % pdb_id)  # 先清空对接文件夹
        os.system(cmdline)
    os.system('cd %s/result/docking_ran/vina && mkdir %s' % (path, pdb_id))  # 创建新的对接文件夹
    ###### 输入文件以及输出路径定义
    
    if os.path.exists(ligand_raw) and os.path.exists(dock_raw) and os.path.exists(protein_raw):
        # 预判蛋白是否准备成功
        cmdline = 'cd %s/result/docking_ran/vina/%s &&' % (path, pdb_id)
        cmdline += 'module purge && '
        cmdline += 'module load mgltools && '
        cmdline += 'prepare_receptor4.py -r %s -o %s.pdbqt -A checkhydrogens -U nphs_lps_waters' % (
            protein_raw, protein_name)
        os.system(cmdline)
        # 预判蛋白是否准备成功
        if os.path.exists('%s/result/docking_ran/vina/%s/%s.pdbqt' % (path, pdb_id, protein_name)):
            cmdline = 'cd %s/result/docking_ran/vina/%s &&' % (path, pdb_id)
            cmdline += 'module purge &&'
            cmdline += 'module load mgltools &&'
            # 准备对接配体
            cmdline += 'prepare_ligand4.py -l %s -o %s.pdbqt -A checkhydrogens ' % (dock_raw, dock_name.split('.')[0])
            os.system(cmdline)
            # 预判配体是否准备成功
            a, b, c, x, y, z = get_mol2_size_center(ligand_raw)  # a, b, c为大小，x,y,z为坐标中心
            if os.path.exists('%s/result/docking_ran/vina/%s/%s.pdbqt' % (path, pdb_id, dock_name.split('.')[0])):
                out_pdbqt = '%s_docked.pdbqt' % dock_name.split('.')[0]
                cmdline = 'cd %s/result/docking_ran/vina/%s &&' % (path, pdb_id)
                cmdline += 'module purge &&'
                cmdline += 'module load vina/1.2.3 &&'
                cmdline += 'vina --ligand %s.pdbqt --receptor %s.pdbqt --center_x %s --center_y %s ' \
                           '--center_z %s --size_x %s --size_y %s --size_z %s --exhaustiveness 32 --energy_range ' \
                           '10000 --num_modes %s --cpu 1 --out %s > %s_docking.log' % (
                               dock_name.split('.')[0], protein_name, x, y, z, a, b, c, num_pose, out_pdbqt, pdb_id)
                os.system(cmdline)
                if os.path.exists(
                        '%s/result/docking_ran/vina/%s/%s_docked.pdbqt' % (path, pdb_id, dock_name.split('.')[0])):
                    # 分割出独立的mol2文件
                    write_mol2(path, pdb_id, dock_name)

                    # 生成*_docking.csv文件
                    # 计算rmsd
                    cmdline = 'cd %s/result/docking_ran/vina/%s &&' % (path, pdb_id)
                    cmdline += 'module purge &&'
                    cmdline += 'module load openbabel &&'
                    cmdline += 'obrms -f %s %s' % (ligand_raw, out_pdbqt.replace('.pdbqt', '.mol2'))  # 根据晶体结构计算rmsd
                    rmsd_res = os.popen(cmdline).readlines()
                    rmsds = [float(_.split()[-1]) for _ in rmsd_res]

                    # 获取energy
                    lines = open('%s/result/docking_ran/vina/%s/%s_docking.log' % (path, pdb_id, pdb_id),
                                 "r").readlines()
                    for idx, line in enumerate(lines):
                        if line.startswith('-----+------------+----------+----------'):
                            st_idx = idx
                    names, energys = [], []
                    for idx, line in enumerate(lines[st_idx + 1:]):
                        idx, energy = float(line.split()[0]), float(line.split()[1])
                        names.append('%s_ligand_vina_%s.mol2' % (pdb_id, int(idx)))
                        energys.append(energy)
                    header = ['name', 'energy', 'rmsd']
                    with open('%s/result/docking_ran/vina/%s/%s_docking.csv' % (path, pdb_id, pdb_id), "w") as f:
                        writer = csv.writer(f)
                        writer.writerow(header)
                        writer.writerows(list(zip(names, energys, rmsds)))


def main(path, pdb_id, num_pose):
    if not os.path.exists('%s/result/docking_ran/vina' % path):
        os.system("cd %s/result/docking_ran && mkdir vina" % path)
    try:
        vina_score(path, pdb_id, num_pose)
    except:
        print('start*************************big error: %s***************************' % pdb_id)


if __name__ == '__main__':
    import time

    st = time.time()
    path = "/home/dejun/workspace/NLDock"
    num_pose = 100
    ligand_files = os.listdir('/home/dejun/workspace/NLDock/docking/ligand/mol2')
    ligand_files = [file for file in ligand_files if file.endswith('_ligand.mol2')]

    pdb_ids = []
    for ligand_file in ligand_files:
        pdb_id = ligand_file[:4]
        if os.path.exists('/home/dejun/workspace/NLDock/result/docking_ran/vina/%s/%s_docking.csv' % (pdb_id, pdb_id)):
            pass
        else:
            # 先清空对接目录
            cmdline = 'rm -rf /home/dejun/workspace/NLDock/result/docking_ran/vina/%s' % pdb_id
            os.system(cmdline)
            pdb_ids.append(pdb_id)

    # pdb_ids = ['1BYJ', '1NEM', '1O9M', '1PBR', '1QD3', '1TOB', '2TOB',
    #            '3SKL', '3SKW', '3SKZ', '3SLM', '3SLQ', '4NYA', '6CC3',
    #            '7ELP', '7ELQ', '7ELR', '7ELS', '7EOG', '7EOK', '7EOL',
    #            '7EOM', '7EON', '7EOO', '7EOP', '7OA3', '7OAV', '7OAW',
    #            '7OAX', '7TZR', '7TZS']
    num_poses = [num_pose for _ in pdb_ids]
    paths = [path for _ in pdb_ids]
    num_cpu = mp.cpu_count()
    pool = mp.Pool(num_cpu)
    pool.starmap_async(main, zip(paths, pdb_ids, num_poses))
    pool.close()
    pool.join()

    # pdb_id = '4NYB'
    # main(path, pdb_id, num_pose)

    end = time.time()
    print('total elapsed time:', end - st, 'S')
