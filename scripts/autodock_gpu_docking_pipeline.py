# ##################################
# autodock/gpu-1.5.3, mgltools/1.5.7
# module purge && source activate /home/dejun/miniconda2/envs/py27 && module load pymol/2.5.4 && python autodock_gpu_docking_pipeline.py
# ========================
import csv
import os
import pymol
from multiprocessing import Pool
from functools import partial
import numpy as np
import time


def get_mol_center_mol2_hf(ligand_mol2):
    x = os.popen(
        "cat %s |sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p |awk '{{print $3}}' | awk '{{x+=$1}} END {{print x/(NR-2)}}'" % ligand_mol2).read()
    y = os.popen(
        "cat %s |sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p |awk '{{print $4}}' | awk '{{y+=$1}} END {{print y/(NR-2)}}'" % ligand_mol2).read()
    z = os.popen(
        "cat %s |sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p |awk '{{print $5}}' | awk '{{z+=$1}} END {{print z/(NR-2)}}'" % ligand_mol2).read()
    print(f'Center of {os.path.basename(ligand_mol2)} is: {x.strip()},{y.strip()},{z.strip()}')
    return float(x.strip()), float(y.strip()), float(z.strip())


def get_mol_size_mol2_hf(ligand_mol2):
    x_max = os.popen(
        "cat %s |sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p |awk '{{print $3}}' |awk 'BEGIN {max = -9999999} {if ($1 > max) max = $1} END {{print max}}'" % ligand_mol2).read()
    x_min = os.popen(
        "cat %s |sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p |awk '{{print $3}}' |awk '{if ($1 != \"\") {print $1}}' |awk 'BEGIN {min = 9999999} {if ($1+0 < min+0) min = $1 fi} END {{print min}}'" % ligand_mol2).read()
    y_max = os.popen(
        "cat %s |sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p |awk '{{print $4}}' |awk 'BEGIN {max = -9999999} {if ($1 > max) max = $1} END {{print max}}'" % ligand_mol2).read()
    y_min = os.popen(
        "cat %s |sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p |awk '{{print $4}}' |awk '{if ($1 != \"\") {print $1}}' |awk 'BEGIN {min = 9999999} {if ($1+0 < min+0) min = $1 fi} END {{print min}}'" % ligand_mol2).read()
    z_max = os.popen(
        "cat %s |sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p |awk '{{print $5}}' |awk 'BEGIN {max = -9999999} {if ($1 > max) max = $1} END {{print max}}'" % ligand_mol2).read()
    z_min = os.popen(
        "cat %s |sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p |awk '{{print $5}}' |awk '{if ($1 != \"\") {print $1}}' |awk 'BEGIN {min = 9999999} {if ($1+0 < min+0) min = $1 fi} END {{print min}}'" % ligand_mol2).read()
    x_size = (float(x_max.strip()) - float(x_min.strip())) * 5
    y_size = (float(y_max.strip()) - float(y_min.strip())) * 5
    z_size = (float(z_max.strip()) - float(z_min.strip())) * 5
    return x_size, y_size, z_size


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
        for idx, line in enumerate(lines[st_idx+1:end_idx]):
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


def write_mol2(i, path, pdb_id, ligand_name):
    cmdline = 'cd %s/result/docking_ran/autodock/%s &&' % (path, pdb_id)
    cmdline += 'module purge &&'
    cmdline += 'module load openbabel &&'
    cmdline += 'obabel -ipdbqt %s_docked_run%s_entity1.pdbqt -omol2 -O %s_autodock_%s.mol2' % (
        ligand_name.split('.')[0], i, ligand_name.split('.')[0], i)
    os.system(cmdline)

    f = open(
        '%s/result/docking_ran/autodock/%s/%s_autodock_%s.mol2' % (
            path, pdb_id, ligand_name.split('.')[0], i), "r")
    lines = f.readlines()
    f.close()
    inx = 0
    content = '%s_autodock_%s.mol2' % (ligand_name.split('.')[0], i)
    for index, line in enumerate(lines):
        if line.startswith("@<TRIPOS>MOLECULE"):
            inx = index + 1
            break
    lines[inx] = content + "\n"
    outfile = open(
        '%s/result/docking_ran/autodock/%s/%s_autodock_%s.mol2' % (
            path, pdb_id, ligand_name.split('.')[0], i),
        "w", newline="")
    for ln in lines:
        outfile.write(ln)
    outfile.close()
    pymol.cmd.load(
        '%s/result/docking_ran/autodock/%s/%s_autodock_%s.mol2' % (
            path, pdb_id, ligand_name.split('.')[0], i))
    pymol.cmd.h_add("all")
    pymol.cmd.save(
        '%s/result/docking_ran/autodock/%s/%s_autodock_%s.mol2' % (
            path, pdb_id, ligand_name.split('.')[0], i))
    pymol.cmd.delete("all")


def autodock_score(path, pdb_id, num_pose):
    ### 输入文件以及输出路径定义
    ligand_name = str(pdb_id + '_ligand.mol2')  # 晶体结构，用于box定义
    ligand_raw = path + '/docking/ligand/mol2/' + ligand_name
    protein_name = pdb_id + '_protein'
    protein_raw = path + '/docking/protein/' + protein_name + '.pdb'

    dock_name = str(pdb_id + '_ligand_random.mol2')
    dock_raw = path + '/docking/ligand_random/mol2/' + dock_name  # 文件位置不一样


    if not os.path.exists((path + '/result/docking_ran/autodock/%s' % pdb_id)):
        os.system('cd %s/result/docking_ran/autodock && mkdir %s' % (path, pdb_id))
    ### 输入文件以及输出路径定义

    if os.path.exists(ligand_raw) and os.path.exists(dock_raw) and os.path.exists(protein_raw):
        cmdline = 'cd %s/result/docking_ran/autodock/%s &&' % (path, pdb_id)
        cmdline += 'module purge && '
        cmdline += 'module load mgltools && '
        cmdline += 'prepare_receptor4.py -r %s -o %s.pdbqt -A checkhydrogens -U nphs_lps_waters' % (
            protein_raw, protein_name)
        os.system(cmdline)
        # 预判蛋白是否准备成功
        if os.path.exists('%s/result/docking_ran/autodock/%s/%s.pdbqt' % (path, pdb_id, protein_name)):
            cmdline = 'cd %s/result/docking_ran/autodock/%s &&' % (path, pdb_id)
            cmdline += 'module purge &&'
            cmdline += 'module load mgltools &&'
            # 准备配体
            cmdline += 'prepare_ligand4.py -l %s -o %s.pdbqt -A checkhydrogens ' % (dock_raw, dock_name.split('.')[0])
            os.system(cmdline)
            # 以晶体构象进行box定义
            a, b, c, x, y, z = get_mol2_size_center(ligand_raw)
            # a, b, c = get_mol_size_mol2_hf(ligand_raw)
            # x, y, z = get_mol_center_mol2_hf(ligand_raw)
            # print('size and center')
            # print(a, b, c, x, y, z)
            # 预判配体是否准备成功
            if os.path.exists('%s/result/docking_ran/autodock/%s/%s.pdbqt' % (path, pdb_id, dock_name.split('.')[0])):
                # 写出sample.gpf文件
                content = 'parameter_file AD4_parameters.dat '

                cmdline = 'cd %s/result/docking_ran/autodock/%s &&' % (path, pdb_id)
                cmdline += 'module purge &&'
                cmdline += 'module load mgltools &&'
                cmdline += "prepare_gpf4.py -l %s.pdbqt -r %s.pdbqt -p npts='%s,%s,%s' -p gridcenter='%s,%s,%s' -p spacing='0.3'" % (
                dock_name.split('.')[0], protein_name, int(a), int(b), int(c), x, y, z)
                os.system(cmdline)
                # 预判gpf是否准备成功
                if os.path.exists('%s/result/docking_ran/autodock/%s/%s.gpf' % (path, pdb_id, protein_name)):
                    cmdline = 'cd %s/result/docking_ran/autodock/%s &&' % (path, pdb_id)
                    cmdline += 'module purge &&'
                    cmdline += 'module load autodock/gpu-1.5.3 &&'
                    cmdline += 'module load cuda/11.2 &&'
                    cmdline += 'autogrid4 -p %s.gpf &&' % protein_name
                    cmdline += 'autodock_gpu_128wi --lfile %s.pdbqt --ffile %s.maps.fld --nrun %s --npdb 1 > %s_docking.log' % (
                        dock_name.split('.')[0],
                        protein_name, num_pose, pdb_id)
                    os.system(cmdline)
                    if os.path.exists('%s/result/docking_ran/autodock/%s/%s.dlg' % (
                            path, pdb_id, dock_name.split('.')[0])) and os.path.exists(
                        '%s/result/docking_ran/autodock/%s/%s_docked_run100_entity1.pdbqt' % (
                        path, pdb_id, dock_name.split('.')[0])):
                        # for 循环改成多进程
                        # for i in range(1, num_pose+1):
                        #     cmdline = 'cd %s/result/docking_ran/autodock/%s &&' % (path, pdb_id)
                        #     cmdline += 'module purge &&'
                        #     cmdline += 'module load openbabel &&'
                        #     cmdline += 'obabel -ipdbqt %s_docked_run%s_entity1.pdbqt -omol2 -O %s_autodock_%s.mol2' % (
                        #     ligand_name.split('.')[0], i, ligand_name.split('.')[0], i)
                        #     os.system(cmdline)
                        #
                        #     f = open(
                        #         '%s/result/docking_ran/autodock/%s/%s_autodock_%s.mol2' % (
                        #         path, pdb_id, ligand_name.split('.')[0], i), "r")
                        #     lines = f.readlines()
                        #     f.close()
                        #     inx = 0
                        #     content = '%s_autodock_%s.mol2' % (ligand_name.split('.')[0], i)
                        #     for index, line in enumerate(lines):
                        #         if line.startswith("@<TRIPOS>MOLECULE"):
                        #             inx = index + 1
                        #             break
                        #     lines[inx] = content + "\n"
                        #     outfile = open(
                        #         '%s/result/docking_ran/autodock/%s/%s_autodock_%s.mol2' % (
                        #         path, pdb_id, ligand_name.split('.')[0], i),
                        #         "w", newline="")
                        #     for ln in lines:
                        #         outfile.write(ln)
                        #     outfile.close()
                        #     pymol.cmd.load(
                        #         '%s/result/docking_ran/autodock/%s/%s_autodock_%s.mol2' % (
                        #         path, pdb_id, ligand_name.split('.')[0], i))
                        #     pymol.cmd.h_add("all")
                        #     pymol.cmd.save(
                        #         '%s/result/docking_ran/autodock/%s/%s_autodock_%s.mol2' % (
                        #         path, pdb_id, ligand_name.split('.')[0], i))
                        #     pymol.cmd.delete("all")

                        # for 循环改成多进程
                        pool = Pool(10)
                        # 写出分子
                        pool.starmap(partial(write_mol2, path=path, pdb_id=pdb_id, ligand_name=dock_name), zip(list(range(1, num_pose+1))))
                        # pool.map(write_mol2, list(range(1, num_pose+1)))
                        pool.close()
                        pool.join()

                        log = open('%s/result/docking_ran/autodock/%s/%s_docking.csv' % (path, pdb_id, pdb_id), "w")
                        mycsv = csv.writer(log)
                        mycsv.writerow(['No.', 'name', 'energy', 'rmsd'])
                        log.close()

                        f = open('%s/result/docking_ran/autodock/%s/%s.dlg' % (path, pdb_id, dock_name.split('.')[0]),
                                 'r')
                        lines = f.readlines()
                        f.close()
                        i = 0
                        for line in lines:
                            if 'RANKING' in line:
                                i += 1
                                numb = line.split()[2]
                                energy = line.split()[3]
                                name = '%s_autodock_%s.mol2' % (dock_name.split('.')[0], numb)
                                cmdline = 'cd %s/result/docking_ran/autodock/%s &&' % (path, pdb_id)
                                cmdline += 'module purge &&'
                                cmdline += 'module load openbabel &&'
                                cmdline += 'obrms %s %s' % (ligand_raw, name)
                                a = os.popen(cmdline).read()
                                rmsd = a.split()[-1]
                                log = open('%s/result/docking_ran/autodock/%s/%s_docking.csv' % (path, pdb_id, pdb_id), "a")
                                mycsv = csv.writer(log)
                                mycsv.writerow([str(i), name, energy, rmsd])
                                log.close()
                        mycsv = open('%s/result/docking_ran/autodock/autodock_run.csv' % path, 'a')
                        mycsvwriter = csv.writer(mycsv)
                        mycsvwriter.writerow([pdb_id, 'dock finish'])
                        mycsv.close()

                        os.system('cd %s/result/docking_ran/autodock/%s && rm *pdbqt *map* *gpf' % (path, pdb_id))

                    else:
                        mycsv = open('%s/result/docking_ran/autodock/autodock_run.csv' % path, 'a')
                        mycsvwriter = csv.writer(mycsv)
                        mycsvwriter.writerow([pdb_id, 'dock failed'])
                        mycsv.close()
                else:
                    mycsv = open('%s/result/docking_ran/autodock/autodock_run.csv' % path, 'a')
                    mycsvwriter = csv.writer(mycsv)
                    mycsvwriter.writerow([pdb_id, 'gpf preparation failed'])
                    mycsv.close()
            else:
                mycsv = open('%s/result/docking_ran/autodock/autodock_run.csv' % path, 'a')
                mycsvwriter = csv.writer(mycsv)
                mycsvwriter.writerow([pdb_id, 'ligand preparation failed'])
                mycsv.close()

        else:
            mycsv = open('%s/result/docking_ran/autodock/autodock_run.csv' % path, 'a')
            mycsvwriter = csv.writer(mycsv)
            mycsvwriter.writerow([pdb_id, 'protein preparation failed'])
            mycsv.close()
    else:
        mycsv = open('%s/result/docking_ran/autodock/autodock_run.csv' % path, 'a')
        mycsvwriter = csv.writer(mycsv)
        mycsvwriter.writerow([pdb_id, 'need one .pdb file and one .mol2 file but got one or none file '])
        mycsv.close()


def main(path, pdb_id, num_pose):
    # 创建CSV
    if not os.path.exists('%s/result/docking_ran/autodock' % path):
        os.system("mkdir -p %s" % '%s/result/docking_ran/autodock' % path)
    if not os.path.exists('%s/result/docking_ran/autodock/autodock_run.csv' % path):
        mycsv = open('%s/result/docking_ran/autodock/autodock_run.csv' % path, 'w')
        mycsvwriter = csv.writer(mycsv)
        mycsvwriter.writerow(
            ['pdbid', 'state']
        )
        mycsv.close()
    try:
        autodock_score(path, pdb_id, num_pose)
    except:
        print('start*************************big error: %s***************************' % pdb_id)


if __name__ == '__main__':
    import time
    st = time.time()
    path = "/home/dejun/workspace/NLDock"
    num_pose = 100

    # 晶体结构
    ligand_files = os.listdir('/home/dejun/workspace/NLDock/docking/ligand/mol2')
    ligand_files = [file for file in ligand_files if file.endswith('_ligand.mol2')]

    for ligand_file in ligand_files:
        pdb_id = ligand_file[:4]
        if os.path.exists('/home/dejun/workspace/NLDock/result/docking_ran/autodock/%s/%s_docking.csv' % (pdb_id, pdb_id)):
            pass
        else:
            # 先清空对接目录
            cmdline = 'rm -rf /home/dejun/workspace/NLDock/result/docking_ran/autodock/%s' % pdb_id
            os.system(cmdline)
            print('start*************************%s***************************' % pdb_id)
            main(path, pdb_id, num_pose)
            print('end*************************%s***************************' % pdb_id)


    # pdbids = ['1BYJ', '1NEM', '1O9M', '1PBR', '1QD3', '1TOB', '2TOB',
    #           '3SKL', '3SKW', '3SKZ', '3SLM', '3SLQ', '4NYA', '6CC3',
    #           '7ELP', '7ELQ', '7ELR', '7ELS', '7EOG', '7EOK', '7EOL',
    #           '7EOM', '7EON', '7EOO', '7EOP', '7OA3', '7OAV', '7OAW',
    #           '7OAX', '7TZR', '7TZS']

    # pdb_id = '1AM0'
    # main(path, pdb_id, num_pose)

    # for pdb_id in pdbids:
    #     # 先清空原来的对接目录
    #     docking_path = '/home/dejun/workspace/NLDock/result/docking_ran/autodock/%s' % pdb_id
    #     if os.path.exists(docking_path):
    #         cmdline = 'rm -rf %s' % docking_path
    #         os.system(cmdline)
    #     main(path, pdb_id, num_pose)
    
    end = time.time()
    print('total elapsed time:', end - st, 'S')
