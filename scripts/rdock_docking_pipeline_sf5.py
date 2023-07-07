##########
### 使用rDock的sf5打分函数对核酸-小分子体系进行对接
##########
import csv
import os
import pandas as pd
import multiprocessing as mp

#  import pymol
num_pose = 100  # 生成构象的数量


def write_file(output_file, outline):
    buffer = open(output_file, 'w')
    buffer.write(outline)
    buffer.close()


# 文件存在，且rmsd不为inf, 最终生成构象数量至少为5个
def valid_docking_check(csv_path):
    flag = True
    if os.path.exists(csv_path):
        csv_data = pd.read_csv(csv_path)
        if len(csv_data) < 5:  # 有可能记录为空
            flag = False
        if all(str(item) == 'inf' for item in csv_data.iloc[:, -1].values):
            flag = False
    else:
        flag = False
    return flag


def rdock_score(pdb_id):
    try:
        print('start*************************%s***************************' % pdb_id)
        ########## 输入文件定义
        ligand_name = str(pdb_id + '_ligand.mol2')
        ligand_raw = path + '/docking/ligand/mol2/' + ligand_name
        protein_name = str(pdb_id + '_protein')
        protein_raw = os.path.join('%s/result/docking_ran/rdock_sf5/%s' % (path, pdb_id),
                                   protein_name + '.mol2')  # 将蛋白文件直接放于对接目录中

        dock_name = str(pdb_id + '_ligand_random.mol2')  # 对接文件定义
        dock_raw = path + '/docking/ligand_random/mol2/' + dock_name
        ########## 输入文件定义

        if os.path.exists(path + '/result/docking_ran/rdock_sf5/%s' % pdb_id):  # 先清空对接目录
            os.system('cd %s/result/docking_ran/rdock_sf5 && rm -rf %s' % (path, pdb_id))
        os.system('cd %s/result/docking_ran/rdock_sf5 && mkdir %s' % (path, pdb_id))  # 创建对接目录

        # 将pdb蛋白结构转换成mol2格式
        if not os.path.exists(protein_raw):
            cmdline = 'module purge && module load schrodinger/2021-1 && structconvert %s %s ' % (
            path + '/docking/protein/' + protein_name + '.pdb', protein_raw)
            os.system(cmdline)
        if os.path.exists(ligand_raw) and os.path.exists(dock_raw) and os.path.exists(protein_raw):
            if not os.path.exists("%s/result/docking_ran/rdock_sf5/%s/%s.sdf" % (path, pdb_id, ligand_name.split(".")[0])):
                cmdline = 'cd %s/result/docking_ran/rdock_sf5/%s &&' % (path, pdb_id)
                cmdline += 'module load openbabel &&'
                cmdline += 'obabel -imol2 %s -osdf -O ./%s.sdf &&' % (ligand_raw, ligand_name.split(".")[0])
                cmdline += 'obabel -imol2 %s -osdf -O ./%s.sdf' % (dock_raw, dock_name.split(".")[0])
                os.system(cmdline)
        if os.path.exists("%s/result/docking_ran/rdock_sf5/%s/%s.sdf" % (path, pdb_id, ligand_name.split(".")[0])):
            if not os.path.exists("%s/result/docking_ran/rdock_sf5/%s/cavity_cav1.grd" % (path, pdb_id)):
                outline = '''RBT_PARAMETER_FILE_V1.00
TITLE %s_dock

RECEPTOR_FILE %s
### RECEPTOR_FLEX 3.0

##################################################################
### CAVITY DEFINITION: REFERENCE LIGAND METHOD
##################################################################
SECTION MAPPER
    SITE_MAPPER RbtLigandSiteMapper
    REF_MOL %s/result/docking_ran/rdock_sf5/%s/%s.sdf
        RADIUS 6.0
        SMALL_SPHERE 1.5
##      LARGE_SPHERE 4.0
        MAX_CAVITIES 1
        MIN_VOLUME 100
        VOL_INCR 0.0
        GRIDSTEP 0.5
END_SECTION

#################################
#CAVITY RESTRAINT PENALTY
#################################
SECTION CAVITY
    SCORING_FUNCTION RbtCavityGridSF
    WEIGHT 1.0
END_SECTION
''' % (pdb_id, protein_raw, path, pdb_id, ligand_name.split(".")[0])
                try:
                    write_file('%s/result/docking_ran/rdock_sf5/%s/cavity.prm' % (path, pdb_id), outline)  # 生成口袋
                    cmdline = 'cd %s/result/docking_ran/rdock_sf5/%s &&' % (path, pdb_id)
                    # cmdline += 'module load rDock &&'
                    cmdline += 'rbcavity -was -d -r cavity.prm'
                    os.system(cmdline)
                except:
                    pass

        # 对构象进行对接
        if os.path.exists("%s/result/docking_ran/rdock_sf5/%s/cavity_cav1.grd" % (path, pdb_id)):
            if not os.path.exists("%s/result/docking_ran/rdock_sf5/%s/stand.txt" % (path, pdb_id)):
                try:
                    cmdline = 'cd %s/result/docking_ran/rdock_sf5/%s &&' % (path, pdb_id)
                    cmdline += 'rbdock -i ./%s.sdf -o stand -r cavity.prm -p /home/dejun/packages/rDock_2022_src/data/scripts/dock_solv.prm -n %s &&' % (
                    dock_name.split(".")[0], num_pose)
                    cmdline += 'sdreport -l stand.sd > stand.txt'
                    os.system(cmdline)
                except:
                    pass

        # 对构象进行打分
        if os.path.exists("%s/result/docking_ran/rdock_sf5/%s/cavity_cav1.grd" % (path, pdb_id)):
            if not os.path.exists("%s/result/docking_ran/rdock_sf5/%s/score.txt" % (path, pdb_id)):
                try:
                    cmdline = 'cd %s/result/docking_ran/rdock_sf5/%s &&' % (path, pdb_id)
                    # cmdline += 'module load rDock &&'
                    cmdline += 'rbdock -i ./%s.sdf -o score -r cavity.prm -p /home/dejun/packages/rDock_2022_src/data/scripts/score_solv.prm -n 1 &&' % (
                    dock_name.split(".")[0])
                    cmdline += 'sdreport -l score.sd > score.txt'
                    os.system(cmdline)
                except:
                    pass

        # 整理对接输出
        if os.path.exists("%s/result/docking_ran/rdock_sf5/%s/stand.txt" % (path, pdb_id)) and os.path.exists(
                "%s/result/docking_ran/rdock_sf5/%s/stand.sd" % (path, pdb_id)):
            energy = list()
            rmsd = list()
            lines = open("%s/result/docking_ran/rdock_sf5/%s/stand.txt" % (path, pdb_id), "r").readlines()
            for line in lines:
                if line.startswith("$SCORE eq"):
                    # print(line.split("eq")[1].strip())
                    energy.append(eval(line.split("eq")[1].strip()))

            cmdline = 'cd %s/result/docking_ran/rdock_sf5/%s &&' % (path, pdb_id)
            cmdline += 'module load openbabel &&'
            cmdline += 'obrms -f %s.sdf stand.sd' % ligand_name.split(".")[0]  # 根据晶体结构计算rmsd
            a = os.popen(cmdline).readlines()
            # lines = a.split('\n')
            for line in a:
                if line.startswith("RMSD"):
                    rmsd.append(line.split()[-1])

            log = open('%s/result/docking_ran/rdock_sf5/%s/%s_docking.csv' % (path, pdb_id, pdb_id), "w")
            mycsv = csv.writer(log)
            mycsv.writerow(['No.', 'energy', 'rmsd'])
            log.close()
            with open('%s/result/docking_ran/rdock_sf5/%s/%s_docking.csv' % (path, pdb_id, pdb_id), "a") as log:
                mycsv = csv.writer(log)
                for i in range(len(energy)):
                    mycsv.writerow([str(i), energy[i], rmsd[i]])
                log.close()
            raw_data = pd.read_csv('%s/result/docking_ran/rdock_sf5/%s/%s_docking.csv' % (path, pdb_id, pdb_id),
                                   index_col=0)
            sort_data = raw_data.sort_values(by=['energy'], ascending=True)
            sort_data.to_csv('%s/result/docking_ran/rdock_sf5/%s/%s_docking.csv' % (path, pdb_id, pdb_id))
        if os.path.exists("%s/result/docking_ran/rdock_sf5/%s/score.txt" % (path, pdb_id)) and os.path.exists(
                "%s/result/docking_ran/rdock_sf5/%s/score.sd" % (path, pdb_id)):
            energy = ''
            lines = open("%s/result/docking_ran/rdock_sf5/%s/score.txt" % (path, pdb_id), "r").readlines()
            for line in lines:
                if line.startswith("$SCORE eq"):
                    # print(line.split("eq")[1].strip())
                    energy = (eval(line.split("eq")[1].strip()))
            with open('%s/result/docking_ran/rdock_sf5/scoring.csv' % path, "a") as log:
                mycsv = csv.writer(log)
                mycsv.writerow([pdb_id, energy])
                log.close()
        print('end*************************%s***************************' % pdb_id)

        # 对生成的decoy构象进行重打分
        # decoys = os.listdir('/home/huifeng/remote/Cyclic_Peptide/work_path/docking_core_ligand/%s' % pdb_id)
        # mycsv = open('%s/result/docking_ran/rdock_sf5/%s/%s_rdock_decoy.csv' % (path, pdb_id, pdb_id), 'w')
        # mycsvwriter = csv.writer(mycsv)
        # mycsvwriter.writerow(
        #     ['ligand_id', 'energy', 'rmsd']
        # )
        # mycsv.close()
        # if not os.path.exists(path + '/result/docking_ran/rdock_sf5/%s/decoy' % pdb_id):
        #     os.system('cd %s/result/docking_ran/rdock_sf5/%s && mkdir decoy' % (path, pdb_id))
        # for ligand in decoys:
        #     if not os.path.exists(path + '/result/docking_ran/rdock_sf5/%s/decoy/%s.sd' % (pdb_id, ligand.split(".")[0])):
        #         try:
        #             cmdline = 'cd %s/result/docking_ran/rdock_sf5/%s &&' % (path, pdb_id)
        #             cmdline += 'module load openbabel &&'
        #             cmdline += 'obabel -imol2 /home/huifeng/remote/Cyclic_Peptide/work_path/docking_core_ligand/%s/%s -osdf -O ./decoy/%s.sdf ' %(pdb_id, ligand, ligand.split(".")[0])
        #             os.system(cmdline)
        #             cmdline = 'cd %s/result/docking_ran/rdock_sf5/%s &&' % (path, pdb_id)
        #             cmdline += 'module load rDock &&'
        #             cmdline += 'rbdock -i decoy/%s.sdf -o decoy/%s -r cavity.prm -p /home/dejun/packages/rDock_2022_src/data/scripts/score_solv.prm -n 1 &&' % (ligand.split(".")[0], ligand.split(".")[0] )
        #             cmdline += 'sdreport -l decoy/%s.sd > decoy/%s.txt' % (ligand.split(".")[0],ligand.split(".")[0])
        #             os.system(cmdline)
        #
        #         except:
        #             pass
        #
        #     if os.path.exists(path + '/result/docking_ran/rdock_sf5/%s/decoy/%s.sd' % (pdb_id, ligand.split(".")[0])) and os.path.exists(path + '/result/docking_ran/rdock_sf5/%s/decoy/%s.txt' % (pdb_id, ligand.split(".")[0])):
        #         energy = ''
        #         lines = open("%s/result/docking_ran/rdock_sf5/%s/decoy/%s.txt"%(path,pdb_id,ligand.split(".")[0]),"r").readlines()
        #         for line in lines:
        #             if line.startswith("$SCORE eq"):
        #                 # print(line.split("eq")[1].strip())
        #                 energy = (eval(line.split("eq")[1].strip()))
        #         rmsd = 0
        #         file_rmsd = open(
        #                 '/home/huifeng/remote/Cyclic_Peptide/work_path/result/docking_randecoy/decoy_refine/%s_rmsd_cluster.csv' % pdb_id,
        #                 'r')
        #         lns = file_rmsd.readlines()
        #         file_rmsd.close()
        #         for ln in lns:
        #             if ln.startswith(ligand):
        #                 rmsd = ln.strip().split(',')[1]
        #
        #         with open('%s/result/docking_ran/rdock_sf5/%s/%s_rdock_decoy.csv' % (path, pdb_id, pdb_id),"a") as log:
        #             mycsv = csv.writer(log)
        #             mycsv.writerow([ligand.split(".")[0], energy, rmsd])
        #             log.close()
    except:
        print('start*************************big error: %s***************************' % pdb_id)


def main(path):
    # 创建CSV
    if not os.path.exists('%s/result/docking_ran/rdock_sf5' % path):
        os.system("mkdir -p %s" % ('%s/result/docking_ran/rdock_sf5' % path))
    if not os.path.exists('%s/result/docking_ran/rdock_sf5/rdock_run.csv' % path):
        mycsv = open('%s/result/docking_ran/rdock_sf5/rdock_run.csv' % path, 'w')
        mycsvwriter = csv.writer(mycsv)
        mycsvwriter.writerow(
            ['pdbid', 'state']
        )
        mycsv.close()
    if not os.path.exists('%s/result/docking_ran/rdock_sf5/rdock_top1.csv' % path):
        for i in ['1', '2', '3']:
            mycsv = open('%s/result/docking_ran/rdock_sf5/rdock_top%s.csv' % (path, i), 'w')
            mycsvwriter = csv.writer(mycsv)
            mycsvwriter.writerow(
                ['pdbid', 'successful_count']
            )
            mycsv.close()
    if not os.path.exists('%s/result/docking_ran/rdock_sf5/scoring.csv' % path):
        log = open('%s/result/docking_ran/rdock_sf5/scoring.csv' % path, "w")
        mycsv = csv.writer(log)
        mycsv.writerow(['pdb_id', 'energy'])
        log.close()

    # 多进程
    protein_names = os.listdir(path + '/docking/protein/')
    protein_names = [name for name in protein_names if name.endswith('_protein.pdb')]
    pdbids = [protein_names[i].split('_')[0] for i in range(len(protein_names))]
    pdbids_ = []
    for pdbid in pdbids:
        # if os.path.exists('%s/result/docking_ran/rdock_sf5/%s/%s_docking.csv' % (path, pdbid, pdbid)):  # 对接成功的蛋白
        if valid_docking_check('%s/result/docking_ran/rdock_sf5/%s/%s_docking.csv' % (path, pdbid, pdbid)):  # 对接成功的蛋白
            pass
        else:  # 没有成功对接的蛋白
            pdbids_.append(pdbid)
            # outdir = '%s/result/docking_ran/rdock_sf5/%s' % (path, pdbid)  # 删除结果目录
            # cmdline = 'rm -rf %s' % outdir
            # os.system(cmdline)
    print(len(pdbids_))
    print(pdbids_)

    # pdbids_ = ['1BYJ', '1NEM', '1O9M', '1PBR', '1QD3', '1TOB', '2TOB',
    #            '3SKL', '3SKW', '3SKZ', '3SLM', '3SLQ', '4NYA', '6CC3',
    #            '7ELP', '7ELQ', '7ELR', '7ELS', '7EOG', '7EOK', '7EOL',
    #            '7EOM', '7EON', '7EOO', '7EOP', '7OA3', '7OAV', '7OAW',
    #            '7OAX', '7TZR', '7TZS']

    pool_num = mp.cpu_count()
    pool = mp.Pool(pool_num)
    pool.starmap_async(rdock_score, zip(pdbids_))
    pool.close()
    pool.join()

    # rdock_score('5SWD')


if __name__ == '__main__':
    import time

    st = time.time()
    path = "/home/dejun/workspace/NLDock/"
    import argparse

    # parser = argparse.ArgumentParser("rdock_docking.")
    # parser.add_argument("pdb_id", action="store", help="input the receptor. E.g:1b6j_C")
    # args = parser.parse_args()
    # pdb_id = args.pdb_id
    main(path)
    end = time.time()
    print('total elapsed time:', end - st, 'S')
