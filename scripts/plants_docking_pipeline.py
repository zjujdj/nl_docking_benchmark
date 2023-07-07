# 利用plants进行对接
import pandas as pd
import numpy as np
import os
import sys
import multiprocessing as mp
import argparse


#########
# module load openbabel/3.1.1
# source activate /home/dejun/conda_env/rtmscore
# python plants_docking_pipeline.py > ./log/plants_docking_pipeline.log
###########


def ligand_prep(ligand_file, target_path):
    cmdline = '%sSPORES_64bit --mode complete %s %s' % (
    plants_path, ligand_file, ligand_file.replace('.mol2', '_sp.mol2'))
    os.system(cmdline)
    cmdline = 'mv %s %s' % (ligand_file.replace('.mol2', '_sp.mol2'), target_path)
    os.system(cmdline)


def protein_prep(mol2_file, target_path):
    # cmdline = 'obabel -ipdb %s -omol2 -O %s' % (pdb_file, pdb_file.replace('.pdb', '.mol2'))  # obabel有时候会自动推断键序，但有时候也不一定
    # cmdline = 'module purge && module load schrodinger/2021-1 && structconvert %s %s' % (pdb_file, pdb_file.replace('.pdb', '.mol2'))
    # os.system(cmdline)
    # 准备蛋白结构 (利用plants自带工具, SPORES_64bit)
    # 利用plants准备好的蛋白文件格式, /%s/%s/pdbid_protein_sp.mol2
    cmdline = '%sSPORES_64bit --mode complete %s %s' % (plants_path, mol2_file, mol2_file.replace('.mol2', '_sp.mol2'))
    os.system(cmdline)
    cmdline = 'mv %s %s' % (mol2_file.replace('.mol2', '_sp.mol2'), target_path)
    os.system(cmdline)


def get_center_from_mol2(ligand_mol2):
    # 根据共晶配体确定坐标
    x = os.popen(
        "cat %s | sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p | awk '{print $3}' | awk '{x+=$1} END {print x/(NR-2)}'" % ligand_mol2).read()
    y = os.popen(
        "cat %s | sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p | awk '{print $4}' | awk '{y+=$1} END {print y/(NR-2)}'" % ligand_mol2).read()
    z = os.popen(
        "cat %s | sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p | awk '{print $5}' | awk '{z+=$1} END {print z/(NR-2)}'" % ligand_mol2).read()
    xyz = [float(x.strip()), float(y.strip()), float(z.strip())]
    return xyz


def get_mol_size(lignd_path):
    x_max = os.popen(
        "cat %s |sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p |awk '{{print $3}}' |awk 'BEGIN {max = 0} {if ($1 > max) max = $1} END {{print max}}'" % lignd_path).read()
    x_min = os.popen(
        "cat %s |sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p |awk '{{print $3}}' |awk '{if ($1 != \"\") {print $1}}' |awk 'BEGIN {min = 9999999} {if ($1+0 < min+0) min = $1 fi} END {{print min}}'" % lignd_path).read()
    x_size = int(float(x_max.strip()) - float(x_min.strip()))
    y_max = os.popen(
        "cat %s |sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p |awk '{{print $4}}' |awk 'BEGIN {max = 0} {if ($1 > max) max = $1} END {{print max}}'" % lignd_path).read()
    y_min = os.popen(
        "cat %s |sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p |awk '{{print $4}}' |awk '{if ($1 != \"\") {print $1}}' |awk 'BEGIN {min = 9999999} {if ($1+0 < min+0) min = $1 fi} END {{print min}}'" % lignd_path).read()
    y_size = int(float(y_max.strip()) - float(y_min.strip()))
    z_max = os.popen(
        "cat %s |sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p |awk '{{print $5}}' |awk 'BEGIN {max = 0} {if ($1 > max) max = $1} END {{print max}}'" % lignd_path).read()
    z_min = os.popen(
        "cat %s |sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p |awk '{{print $5}}' |awk '{if ($1 != \"\") {print $1}}' |awk 'BEGIN {min = 9999999} {if ($1+0 < min+0) min = $1 fi} END {{print min}}'" % lignd_path).read()
    z_size = int(float(z_max.strip()) - float(z_min.strip()))
    xyz = [x_size, y_size, z_size]
    return xyz


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
    x_size = (x_max - x_min)
    y_size = (y_max - y_min)
    z_size = (z_max - z_min)
    x_center, y_center, z_center = sum(x_coord) / n_atoms, sum(y_coord) / n_atoms, sum(z_coord) / n_atoms
    return x_size, y_size, z_size, x_center, y_center, z_center


def write_config_file(protein_file, ligand_file, config_file, output_dir, x, y, z, x_size, y_size, z_size):
    content = f'''# scoring function and search settings
# plp|plp95|chemplp
scoring_function chemplp
search_speed speed1
# input
protein_file {protein_file}
ligand_file {ligand_file}
# output
output_dir {output_dir}
# write single mol2 files (e.g. for RMSD calculation)
write_multi_mol2 1
write_per_atom_scores 0
write_protein_conformations 0
write_protein_bindingsite 0
# binding site definition
bindingsite_center {x} {y} {z}
bindingsite_radius {max([x_size, y_size, z_size]) / 1.5}
# cluster algorithm
cluster_structures {num_pose}
cluster_rmsd 2.0
'''
    with open(config_file, 'w') as f:
        f.write(content)


def plants_docking(config_file):
    cmd = f'timeout 20m %sPLANTS1.2_64bit --mode screen %s' % (plants_path, config_file)
    # cmd = f'%sPLANTS1.2_64bit --mode screen %s' % (plants_path, config_file)
    os.system(cmd)


def plants_docking_pipeline(pdb_id):
    try:
        print('start*************************%s***************************' % pdb_id)
        ######### 输入文件定义
        protein_file = '/home/dejun/workspace/NLDock/docking/protein/%s_protein.mol2' % pdb_id
        ligand_file = '/home/dejun/workspace/NLDock/docking/ligand/mol2/%s_ligand.mol2' % pdb_id

        dock_file = '/home/dejun/workspace/NLDock/docking/ligand_random/mol2/%s_ligand_random.mol2' % pdb_id
        ######### 输入文件定义

        res_path = '/home/dejun/workspace/NLDock/result/docking_ran/plants'
        if not os.path.exists(res_path):
            cmdline = 'mkdir -p %s' % res_path
            os.system(cmdline)

        if os.path.exists('%s/%s' % (res_path, pdb_id)):  # 先清空对接目录
            cmdline = 'rm -rf %s' % '%s/%s' % (res_path, pdb_id)
            os.system(cmdline)
        cmdline = 'cd %s && mkdir %s' % (res_path, pdb_id)  # 新建对接目录
        os.system(cmdline)

        # config file写出路径
        output_dir = '%s/%s/docking' % (res_path, pdb_id)
        # 准备对接配体
        ligand_prep(dock_file, '%s/%s' % (res_path, pdb_id))
        # 准备蛋白
        protein_prep(protein_file, '%s/%s' % (res_path, pdb_id))
        a, b, c, x, y, z = get_mol2_size_center(ligand_file)  # 根据对接晶体结构获取对接box信息
        print('%s: size and center' % pdb_id)
        print(a, b, c, x, y, z)
        # print(xyz, xyz_size)
        config_file = '%s/%s/%s.conf' % (res_path, pdb_id, pdb_id)
        # 写出对接的config文件
        write_config_file(
            os.path.join('%s/%s' % (res_path, pdb_id), os.path.basename(protein_file.replace('.mol2', '_sp.mol2'))),
            os.path.join('%s/%s' % (res_path, pdb_id),
                         os.path.basename(dock_file.replace('.mol2', '_sp.mol2')))
            , config_file, output_dir, x, y, z, a, b, c)
        # 分子对接
        plants_docking(config_file)

        if os.path.exists('%s/docked_ligands.mol2' % output_dir):
            cmdline = 'module load openbabel/3.1.1 && obrms -f %s %s' % (ligand_file, '%s/docked_ligands.mol2' % output_dir)  # 根据晶体结构计算rmsd
            rmsd_res = os.popen(cmdline).read()
            rmsd_res = rmsd_res.split('\n')[:-1]  # 根据'\n'分割出的list最后一个元素为''
            rmsds = [eval(line.split()[-1]) for line in rmsd_res]

            csv_file = pd.read_csv('%s/ranking.csv' % output_dir)
            csv_file['rmsd'] = rmsds

            def get_name(title):
                id = int(title.split()[-1].split('_')[-1])
                name = '%s_ligand_random_plants_%s.mol2' % (pdb_id, id)  ################################################
                # name = '%s_ligand_plants_%s.mol2' % (pdb_id, id)
                return name

            csv_file['name'] = csv_file['LIGAND_ENTRY'].apply(get_name)
            csv_file['energy'] = csv_file['TOTAL_SCORE']
            csv_file[['name', 'energy', 'rmsd']].to_csv('%s/%s_docking.csv' % (output_dir, pdb_id), index=False)
            print('end*************************%s***************************' % pdb_id)
    except:
        print('*************************big error: %s***************************' % pdb_id)


if __name__ == '__main__':
    import time

    st = time.time()
    #
    plants_path = '/home/dejun/workspace/metal_binding/ign_metal_prediction_pipeline/plants/'
    num_pose = 100
    path = "/home/dejun/workspace/NLDock"

    protein_names = os.listdir(path + '/docking/protein/')
    protein_names = [name for name in protein_names if name.endswith('_protein.pdb')]
    pdbids = [protein_names[i].split('_')[0] for i in range(len(protein_names))]
    pdbids_ = []
    for pdbid in pdbids:
        if os.path.exists('%s/result/docking_ran/plants/%s/docking/%s_docking.csv' % (path, pdbid, pdbid)):  # 对接成功的蛋白
            pass
        else:  # 没有成功对接的蛋白
            pdbids_.append(pdbid)
            outdir = '%s/result/docking_ran/plants/%s' % (path, pdbid)  # 删除结果目录
            cmdline = 'rm -rf %s' % outdir
            os.system(cmdline)

    # pdbids_ = ['1BYJ', '1NEM', '1O9M', '1PBR', '1QD3', '1TOB', '2TOB',
    #            '3SKL', '3SKW', '3SKZ', '3SLM', '3SLQ', '4NYA', '6CC3',
    #            '7ELP', '7ELQ', '7ELR', '7ELS', '7EOG', '7EOK', '7EOL',
    #            '7EOM', '7EON', '7EOO', '7EOP', '7OA3', '7OAV', '7OAW',
    #            '7OAX', '7TZR', '7TZS']

    pool_num = mp.cpu_count()
    pool = mp.Pool(pool_num)
    pool.starmap_async(plants_docking_pipeline, zip(pdbids_))
    pool.close()
    pool.join()

    # # 测试单个蛋白
    # pdbid = '7TZT'
    # outdir = '%s/result/docking_ran/plants/%s' % (path, pdbid)  # 删除结果目录
    # cmdline = 'rm -rf %s' % outdir
    # os.system(cmdline)
    # plants_docking_pipeline('7TZT')

    end = time.time()
    print('total elapsed time:', end - st, 'S')
