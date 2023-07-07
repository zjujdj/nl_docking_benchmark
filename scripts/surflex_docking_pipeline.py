# tripos/sybyl/x2.1.1
import csv
import os
import time
from multiprocessing.dummy import Pool
import multiprocessing as mp
num_pose = 100  # 最多产生多少个构象


def get_pocket(reDOCKINGid):
    ligand_name = reDOCKINGid + '_ligand.mol2'
    protein_name = reDOCKINGid + '_protein'

    cmdline = 'cd %s/result/docking_ran/surflex/%s &&' % (path, reDOCKINGid)
    cmdline += 'module load tripos/sybyl/x2.1.1 &&'
    cmdline += 'surflex-dock.exe proto %s/docking/ligand/mol2/%s ./%s.mol2 p1 && ' % (path, ligand_name, protein_name)
    cmdline += 'rm -rf p1-*.pdb'
    os.system(cmdline)


def write_file(output_filename, outline):
    buffer = open(output_filename, 'w')
    buffer.write(outline)
    buffer.close()


def surflex_decoy(reDOCKINGid):
    try:
        print('start*************************%s***************************' % reDOCKINGid)
        ########## 输入文件定义
        ligand_name = reDOCKINGid + '_ligand.mol2'
        protein_name = reDOCKINGid + '_protein'
        dock_name = reDOCKINGid + '_ligand_random.mol2'
        ########## 输入文件定义


        # 清空对接目录
        cmdline = 'cd %s/result/docking_ran/surflex && ' % path
        cmdline += 'rm -rf %s' % reDOCKINGid
        os.system(cmdline)

        # 创建对接目录
        cmdline = 'cd %s/result/docking_ran/surflex && ' % path
        cmdline += 'mkdir %s' % reDOCKINGid
        os.system(cmdline)

        # 转换蛋白格式, pdb2mol2
        cmdline = 'cd %s/docking/protein &&' % path
        cmdline += 'module load openeye &&'
        cmdline += 'convert.py ./%s.pdb %s/result/docking_ran/surflex/%s/%s.mol2 ' % (protein_name, path, reDOCKINGid, protein_name)
        os.system(cmdline)

        if os.path.exists('%s/result/docking_ran/surflex/%s/%s.mol2' % (path, reDOCKINGid, protein_name)):
            # 获取蛋白质口袋信息
            get_pocket(reDOCKINGid)
            if os.path.exists('%s/result/docking_ran/surflex/%s/p1-protomol.mol2' % (path, reDOCKINGid)):
                write_file('%s/result/docking_ran/surflex/%s/list' % (path, reDOCKINGid), '%s/docking/ligand_random/mol2/%s' % (path, dock_name))  ############### 需要对接分子的list
                try:
                    # cmdline = 'cd %s/result/docking_ran/surflex/%s && ' % (path, reDOCKINGid)
                    # cmdline += 'module load tripos/sybyl/x2.1.1 &&'
                    # cmdline += 'surflex-dock.exe -ndock_final 100 -multistart 2 -maxrot 100 -protodthresh 15.0 +fastsearch -maxconfs 20 dock_list ./list ./p1-protomol.mol2 ./%s.mol2 %s_2_loggeom' % (protein_name,reDOCKINGid)
                    # os.system(cmdline)
                    # cmdline = 'cd %s/result/docking_ran/surflex/%s && '% (path, reDOCKINGid)
                    # cmdline += 'module load tripos/sybyl/x2.1.1 &&'
                    # cmdline += 'surflex-dock.exe -ndock_final 100 -multistart 5 -maxrot 100 -protodthresh 15.0 +fastsearch -maxconfs 20 dock_list ./list ./p1-protomol.mol2 ./%s.mol2 %s_5_loggeom' % (protein_name,reDOCKINGid)
                    # os.system(cmdline)
                    cmdline = 'cd %s/result/docking_ran/surflex/%s && ' % (path, reDOCKINGid)
                    cmdline += 'module load tripos/sybyl/x2.1.1 &&'
                    cmdline += 'timeout 20m surflex-dock.exe -ndock_final %s -multistart 10 -maxrot 100 -protodthresh 15.0 +fastsearch -maxconfs 20 dock_list ./list ./p1-protomol.mol2 ./%s.mol2 %s_10_loggeom' % (num_pose, protein_name, reDOCKINGid)
                    os.system(cmdline)
                except:
                    # mycsv = open('%s/result/docking_ran/surflex/surflex_decoy_2_out.csv' % path, 'a')
                    # mycsvwriter = csv.writer(mycsv)
                    # mycsvwriter.writerow([reDOCKINGid, "decoy_genate_failed"])
                    # mycsv.close()
                    # mycsv = open('%s/result/docking_ran/surflex/surflex_decoy_5_out.csv' % path, 'a')
                    # mycsvwriter = csv.writer(mycsv)
                    # mycsvwriter.writerow([reDOCKINGid, "decoy_genate_failed"])
                    # mycsv.close()
                    mycsv = open('%s/result/docking_ran/surflex/surflex_decoy_10_out.csv' % path, 'a')
                    mycsvwriter = csv.writer(mycsv)
                    mycsvwriter.writerow([reDOCKINGid, "decoy_genate_failed"])
                    mycsv.close()

                # 计算RMSD
                gen_pose_file = '%s/result/docking_ran/surflex/%s/%s_10_loggeom-results.mol2' % (path, reDOCKINGid, reDOCKINGid)
                lig_pose_file = '%s/docking/ligand_random/mol2/%s_ligand.mol2' % (path, reDOCKINGid)
                rmsd_out_file = '%s/result/docking_ran/surflex/%s/%s_rmsd.txt' % (path, reDOCKINGid, reDOCKINGid)
                if os.path.exists(gen_pose_file):
                    cmdline = 'obrms -f %s %s > %s' % (lig_pose_file, gen_pose_file, rmsd_out_file)
                    os.system(cmdline)

                    # 提取RMSD和打分值到csv文件中
                    if not os.path.exists('%s/result/docking_ran/surflex/%s/%s_docking.csv' % (path, reDOCKINGid, reDOCKINGid)):
                        mycsv = open('%s/result/docking_ran/surflex/%s/%s_docking.csv' % (path, reDOCKINGid, reDOCKINGid), 'w')
                        mycsvwriter = csv.writer(mycsv)
                        mycsvwriter.writerow(['name', 'energy', 'rmsd'])
                        mycsv.close()

                        lns = open(rmsd_out_file, 'r').readlines()
                        file = open('%s/result/docking_ran/surflex/%s/%s_10_loggeom' % (path, reDOCKINGid, reDOCKINGid), 'r')
                        line = file.readline()
                        i = 0  # 第0行不是记录的打分值
                        while True:
                            line = file.readline()
                            if not line:
                                break
                            energy = line.strip().split()[1]
                            # energy = eval(line.split()[3])
                            i += 1
                            rmsd = 'inf'
                            # lns = open(rmsd_out_file, 'r').readlines()
                            # for ln in lns:
                            #     if 'surflex_10_%s' % (str(i - 1).zfill(3)) in ln:
                            #         rmsd = ln.strip().split()[-1]
                            rmsd = eval(lns[i-1].strip().split()[-1])
                            mycsv = open('%s/result/docking_ran/surflex/%s/%s_docking.csv' % (path, reDOCKINGid, reDOCKINGid), 'a')
                            mycsvwriter = csv.writer(mycsv)
                            mycsvwriter.writerow(['%s_surflex_%s.mol2' % (reDOCKINGid, str(i - 1).zfill(3)), energy, rmsd])
                            mycsv.close()
                        file.close()

                # # # 用不同的子目录保存对接结果
                # cmdline = "cd %s/result/docking_ran/surflex/%s && " % (path, reDOCKINGid)
                # cmdline += "mv %s_*_loggeom-results.mol2 %s_*_loggeom -t ../ &&" %(reDOCKINGid, reDOCKINGid)
                # cmdline += "cd ../ &&"
                # cmdline += "rm -r %s" %(reDOCKINGid)
                # os.system(cmdline)


            else:
                # mycsv = open('%s/result/docking_ran/surflex/surflex_decoy_2_out.csv' % path, 'a')
                # mycsvwriter = csv.writer(mycsv)
                # mycsvwriter.writerow([reDOCKINGid, "no_pocket"])
                # mycsv.close()
                # mycsv = open('%s/result/docking_ran/surflex/surflex_decoy_5_out.csv' % path, 'a')
                # mycsvwriter = csv.writer(mycsv)
                # mycsvwriter.writerow([reDOCKINGid, "no_pocket"])
                # mycsv.close()
                mycsv = open('%s/result/docking_ran/surflex/surflex_decoy_10_out.csv' % path, 'a')
                mycsvwriter = csv.writer(mycsv)
                mycsvwriter.writerow([reDOCKINGid, "no_pocket"])
                mycsv.close()

        else:
            # mycsv = open('%s/result/docking_ran/surflex/surflex_decoy_2_out.csv' % path, 'a')
            # mycsvwriter = csv.writer(mycsv)
            # mycsvwriter.writerow([reDOCKINGid,"protein_convert_failed"])
            # mycsv.close()
            # mycsv = open('%s/result/docking_ran/surflex/surflex_decoy_5_out.csv' % path, 'a')
            # mycsvwriter = csv.writer(mycsv)
            # mycsvwriter.writerow([reDOCKINGid,"protein_convert_failed"])
            # mycsv.close()
            mycsv = open('%s/result/docking_ran/surflex/surflex_decoy_10_out.csv' % path, 'a')
            mycsvwriter = csv.writer(mycsv)
            mycsvwriter.writerow([reDOCKINGid, "protein_convert_failed"])
            mycsv.close()
        print('end*************************%s***************************' % reDOCKINGid)
    except:
        print('start*************************big error: %s***************************' % reDOCKINGid)


def main(path):
    # 创建CSV
    if not os.path.exists('%s/result/docking_ran/surflex' % path):
        os.system("mkdir -p %s" % '%s/result/docking_ran/surflex' % path)
    # mycsv = open('%s/result/docking_ran/surflex/surflex_decoy_2_out.csv' % path , 'w')
    # mycsvwriter = csv.writer(mycsv)
    # mycsvwriter.writerow(['name', 'state'])
    # mycsv.close()
    # mycsv = open('%s/result/docking_ran/surflex/surflex_decoy_5_out.csv' % path , 'w')
    # mycsvwriter = csv.writer(mycsv)
    # mycsvwriter.writerow(['name', 'state'])
    # mycsv.close()
    if os.path.exists('%s/result/docking_ran/surflex/surflex_decoy_10_out.csv' % path):
        pass
    else:
        mycsv = open('%s/result/docking_ran/surflex/surflex_decoy_10_out.csv' % path, 'w')
        mycsvwriter = csv.writer(mycsv)
        mycsvwriter.writerow(['name', 'state'])
        mycsv.close()

    # 多进程
    protein_names = os.listdir(path + '/docking/protein/')
    protein_names = [name for name in protein_names if name.endswith('_protein.pdb')]
    pdbids = [protein_names[i].split('_')[0] for i in range(len(protein_names))]
    pdbids_ = []
    for pdbid in pdbids:
        if os.path.exists('%s/result/docking_ran/surflex/%s/%s_docking.csv' % (path, pdbid, pdbid)):  # 对接成功的蛋白
            pass
        else:  # 没有成功对接的蛋白
            pdbids_.append(pdbid)
            outdir = '%s/result/docking_ran/surflex/%s' % (path, pdbid)  # 删除结果目录
            cmdline = 'rm -rf %s' % outdir
            os.system(cmdline)
    # pdbids_ = ['1BYJ', '1NEM', '1O9M', '1PBR', '1QD3', '1TOB', '2TOB',
    #            '3SKL', '3SKW', '3SKZ', '3SLM', '3SLQ', '4NYA', '6CC3',
    #            '7ELP', '7ELQ', '7ELR', '7ELS', '7EOG', '7EOK', '7EOL',
    #            '7EOM', '7EON', '7EOO', '7EOP', '7OA3', '7OAV', '7OAW',
    #            '7OAX', '7TZR', '7TZS']
    pool_num = mp.cpu_count()
    pool = mp.Pool(pool_num)
    pool.starmap_async(surflex_decoy, zip(pdbids_))
    pool.close()
    pool.join()

    # surflex_decoy("5Z1H")


if __name__ == '__main__':
    import time
    st = time.time()
    # path = sys.argv[1]
    path = "/home/dejun/workspace/NLDock"
    pool_num = mp.cpu_count()
    main(path)
    end = time.time()
    print('total elapsed time:', end - st, 'S')

