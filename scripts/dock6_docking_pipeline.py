# /home/dejun/workspace/NLDock/pymol/bin/python dock6_docking_pipeline.py
import os
import re
import sys
from functools import partial
import multiprocessing as mp
import pandas as pd


def read_file(file_name):
    with open(file_name, 'r') as f:
        lines = f.readlines()
    return lines


def write_file(outline, file_name):
    with open(file_name, 'w') as f:
        f.write(outline)


def run_dock6(pdb_id):
    try:
        ############## 输入文件定义
        crystal_ligand = '/home/dejun/workspace/NLDock/docking/ligand/mol2/%s_ligand.mol2' % pdb_id
        # receptor_file_ = '/home/dejun/workspace/NLDock/docking/protein/%s_protein.pdb' % pdb_id
        receptor_file = '/home/dejun/workspace/NLDock/docking/protein/%s_protein.mol2' % pdb_id
        # 先利用structconvert将蛋白pdb转换成mol2
        # cmdline = 'module purge && module load schrodinger/2021-1 && structconvert %s %s' % (receptor_file_, receptor_file_.replace('.pdb', '_rec.mol2'))
        # os.system(cmdline)
        ligand_file = '/home/dejun/workspace/NLDock/docking/ligand_random/mol2/%s_ligand_random.mol2' % pdb_id  # 对接文件
        ############## 输入文件定义
        
        dock6_home = '/opt/dock/6.10'
        docking_path = '/home/dejun/workspace/NLDock/result/docking_ran/dock6/%s' % pdb_id
        if os.path.exists(docking_path):
            cmdline = 'rm -rf  %s' % docking_path  # 先清空对接文件目录
            os.system(cmdline)
            
        cmdline = 'mkdir -p %s' % docking_path   # 创建新的对接目录
        os.system(cmdline)
        # 将转换好的mol2受体文件移动到docking_path
        # cmdline = 'mv %s %s' % (receptor_file_.replace('.pdb', '_rec.mol2'), docking_path)
        # os.system(cmdline)
        # receptor_file = os.path.join(docking_path, os.path.basename(receptor_file_.replace('.pdb', '_rec.mol2')))

        receptor_name = os.path.basename(receptor_file)
        ligand_name = os.path.basename(ligand_file)

        # write_dms.py
        outline = 'from chimera import openModels, Molecule, runCommand, MSMSModel\n'
        outline += 'from WriteDMS import writeDMS\n'
        outline += '\n'
        outline += 'models = openModels.list(modelTypes=[Molecule])\n'
        outline += "runCommand('delete element.H')\n"
        outline += "runCommand('surf')\n"
        outline += 'surf = openModels.list(modelTypes=[MSMSModel])[0]\n'
        outline += "writeDMS(surf, '%s/receptor.dms')\n" % docking_path
        write_file(outline, os.path.join(docking_path, 'write_dms.py'))

        # INSPH
        outline = '%s/receptor.dms\n' % docking_path
        outline += 'R\n'
        outline += 'X\n'
        outline += '0.0\n'
        outline += '4.0\n'
        outline += '1.4\n'
        outline += '%s/receptor.sph\n' % docking_path
        write_file(outline, os.path.join(docking_path, 'INSPH'))

        # box.in
        outline = 'Y\n'
        outline += '8.0\n'
        outline += '%s/selected_spheres.sph\n' % docking_path
        outline += '1\n'
        outline += '%s/box.pdb' % docking_path
        write_file(outline, os.path.join(docking_path, 'box.in'))

        # grid.in
        outline = 'compute_grids            yes\n'
        outline += 'grid_spacing            0.3\n'
        outline += 'output_molecule         no\n'
        outline += 'contact_score           no\n'
        outline += 'energy_score            yes\n'
        outline += 'energy_cutoff_distance  9999\n'
        outline += 'atom_model              a\n'
        outline += 'attractive_exponent     6\n'
        outline += 'repulsive_exponent      12\n'
        outline += 'distance_dielectric     yes\n'
        outline += 'allow_non_integral_charges     yes\n'
        outline += 'dielectric_factor       4\n'
        outline += 'bump_filter             yes\n'
        outline += 'bump_overlap            0.75\n'
        outline += 'receptor_file           %s\n' % receptor_file
        outline += 'box_file                /%s/box.pdb\n' % docking_path
        outline += 'vdw_definition_file     %s/parameters/vdw_AMBER_parm99.defn\n' % dock6_home
        outline += 'score_grid_prefix       grid'
        write_file(outline, os.path.join(docking_path, 'grid.in'))

        # dock.in
        outline = 'ligand_atom_file                       %s\n' % ligand_file
        outline += 'limit_max_ligands                     no\n'
        outline += 'skip_molecule                         no\n'
        outline += 'read_mol_solvation                    no\n'
        outline += 'calculate_rmsd                        no\n'
        outline += 'use_database_filter                   no\n'
        outline += 'orient_ligand                         yes\n'
        outline += 'automated_matching                    yes\n'
        outline += 'receptor_site_file                    %s/selected_spheres.sph\n' % docking_path
        outline += 'max_orientations                      50000\n'
        outline += 'critical_points                       no\n'
        outline += 'chemical_matching                     no\n'
        outline += 'use_ligand_spheres                    no\n'
        outline += 'use_internal_energy                   yes\n'
        outline += 'internal_energy_rep_exp               12\n'
        outline += 'flexible_ligand                       yes\n'
        outline += 'user_specified_anchor                 no\n'
        outline += 'limit_max_anchors                     no\n'
        outline += 'min_anchor_size                       6\n'
        outline += 'pruning_use_clustering                yes\n'
        outline += 'pruning_max_orients                   1000\n'
        outline += 'pruning_clustering_cutoff             100\n'
        outline += 'pruning_conformer_score_cutoff        100.0\n'
        outline += 'use_clash_overlap                     no\n'
        outline += 'write_growth_tree                     no\n'
        outline += 'bump_filter                           no\n'
        outline += 'score_molecules                       yes\n'
        outline += 'contact_score_primary                 no\n'
        outline += 'contact_score_secondary               no\n'
        outline += 'grid_score_primary                    yes\n'
        outline += 'grid_score_secondary                  no\n'
        outline += 'grid_score_rep_rad_scale              1\n'
        outline += 'grid_score_vdw_scale                  1\n'
        outline += 'grid_score_es_scale                   1\n'
        outline += 'grid_score_grid_prefix                grid\n'
        outline += 'multigrid_score_secondary             no\n'
        outline += 'dock3.5_score_secondary               no\n'
        outline += 'continuous_score_secondary            no\n'
        outline += 'descriptor_score_secondary            no\n'
        outline += 'gbsa_zou_score_secondary              no\n'
        outline += 'gbsa_hawkins_score_secondary          no\n'
        outline += 'SASA_descriptor_score_secondary       no\n'
        outline += 'amber_score_secondary                 no\n'
        outline += 'minimize_ligand                       yes\n'
        outline += 'minimize_anchor                       yes\n'
        outline += 'minimize_flexible_growth              yes\n'
        outline += 'use_advanced_simplex_parameters       no\n'
        outline += 'simplex_max_cycles                    1\n'
        outline += 'simplex_score_converge                0.1\n'
        outline += 'simplex_cycle_converge                1.0\n'
        outline += 'simplex_trans_step                    1.0\n'
        outline += 'simplex_rot_step                      0.1\n'
        outline += 'simplex_tors_step                     10.0\n'
        outline += 'simplex_anchor_max_iterations         500\n'
        outline += 'simplex_grow_max_iterations           500\n'
        outline += 'simplex_grow_tors_premin_iterations   0\n'
        outline += 'simplex_random_seed                   0\n'
        outline += 'simplex_restraint_min                 no\n'
        outline += 'atom_model                            all\n'
        outline += 'vdw_defn_file                         %s/parameters/vdw_AMBER_parm99.defn\n' % dock6_home
        outline += 'flex_defn_file                        %s/parameters/flex.defn\n' % dock6_home
        outline += 'flex_drive_file                       %s/parameters/flex_drive.tbl\n' % dock6_home
        outline += 'ligand_outfile_prefix                 docked\n'
        outline += 'write_orientations                    no\n'
        outline += 'num_scored_conformers                 100\n'
        outline += 'write_conformations                   yes\n'
        outline += 'cluster_conformations                 yes\n'
        outline += 'cluster_rmsd_threshold                2.0\n'
        outline += 'rank_ligands                          no'
        write_file(outline, os.path.join(docking_path, 'dock.in'))

        # # run grid and dock
        # cmdline = 'export PATH=%s/bin:${PATH} && ' % dock6_home
        # cmdline += 'unset AMBERHOME && '
        # cmdline += 'chimera --nogui %s write_dms.py > write_dms.log && ' % receptor_file
        # cmdline += 'sphgen_cpp -i receptor.dms -o receptor.sph > outsph && '
        # cmdline += 'sphere_selector receptor.sph %s 10.0 && ' % crystal_ligand
        # cmdline += 'showbox < box.in > showbox.log && '
        # cmdline += 'grid -i grid.in -o grid.out && '
        # cmdline += 'dock6 -i dock.in -o dock.out'

        # run grid and dock
        cmdline = 'cd %s && export PATH=%s/bin:${PATH} && ' % (docking_path, dock6_home)
        cmdline += 'unset AMBERHOME && '
        cmdline += 'module load chimera/1.15 && module load dock/6.10 && chimera --nogui %s %s/write_dms.py > %s/write_dms.log && ' % (receptor_file, docking_path, docking_path)
        cmdline += 'sphgen -i %s/INSPH && ' % docking_path
        cmdline += 'sphere_selector %s/receptor.sph %s 10.0 && ' % (docking_path, crystal_ligand)
        cmdline += 'showbox < %s/box.in > %s/showbox.log && ' % (docking_path, docking_path)
        cmdline += 'grid -i %s/grid.in -o %s/grid.out && ' % (docking_path, docking_path)
        cmdline += 'dock6 -i %s/dock.in -o %s/dock.out' % (docking_path, docking_path)

        os.system(cmdline)

        # 写出单个分子以及每个分子的得分和rmsd
        with open('%s/docked_scored.mol2' % docking_path, 'r') as f:
            mol_lines = f.readlines()

        names, mark_idxs = [], []
        for idx, line in enumerate(mol_lines):
            if line.startswith('##########                                Name:'):
                mark_idxs.append(idx)
        # write molecules
        for idx, mol_idx in enumerate(mark_idxs[:-1]):
            name = '%s_ligand_random_dock6_%s.mol2' % (pdb_id, idx + 1)  ########################################################
            # name = '%s_ligand_dock6_%s.mol2' % (pdb_id, idx + 1)
            names.append(name)
            with open('%s/%s' % (docking_path, name), 'w') as f:
                f.write(''.join(mol_lines[mol_idx:mark_idxs[idx + 1]]))

        # 最后一个分子
        name = '%s_ligand_random_dock6_%s.mol2' % (pdb_id, len(mark_idxs))  ########################################################
        # name = '%s_ligand_dock6_%s.mol2' % (pdb_id, len(mark_idxs))
        names.append(name)
        with open('%s/%s' % (docking_path, name), 'w') as f:
            f.write(''.join(mol_lines[mark_idxs[-1]:]))

        # 记录能量得分
        scores = []
        for line in mol_lines:
            if line.startswith('##########                          Grid_Score:'):
                scores.append(float(line.split()[-1]))
        # rmsd
        cmdline = 'module load openbabel/3.1.1 && obrms -f %s %s/docked_scored.mol2' % (crystal_ligand, docking_path)
        lines = os.popen(cmdline).readlines()
        rmsds = [line.split()[-1] for line in lines]

        csv_file = pd.DataFrame({'name': names, 'energy': scores, 'rmsd': rmsds})
        csv_file.to_csv('%s/%s_docking.csv' % (docking_path, pdb_id), index=False)
    except:
        print('***********error:%s*******************' % pdb_id)


if __name__ == '__main__':
    import time

    st = time.time()
    path = "/home/dejun/workspace/NLDock"
    ligand_files = os.listdir('/home/dejun/workspace/NLDock/docking/ligand/mol2')
    ligand_files = [file for file in ligand_files if file.endswith('_ligand.mol2')]

    pdb_ids = []
    for ligand_file in ligand_files:
        pdb_id = ligand_file[:4]
        if os.path.exists(
                '/home/dejun/workspace/NLDock/result/docking_ran/dock6/%s/%s_docking.csv' % (pdb_id, pdb_id)):
            pass
        else:
            pdb_ids.append(pdb_id)
            # 需要先删除目录
            cmdline = 'rm -rf /home/dejun/workspace/NLDock/result/docking_ran/dock6/%s' % pdb_id
            os.system(cmdline)
    print(len(pdb_ids))
    print(pdb_ids)

    # pdbids = ['1BYJ', '1NEM', '1O9M', '1PBR', '1QD3', '1TOB', '2TOB',
    #           '3SKL', '3SKW', '3SKZ', '3SLM', '3SLQ', '4NYA', '6CC3',
    #           '7ELP', '7ELQ', '7ELR', '7ELS', '7EOG', '7EOK', '7EOL',
    #           '7EOM', '7EON', '7EOO', '7EOP', '7OA3', '7OAV', '7OAW',
    #           '7OAX', '7TZR', '7TZS']
    
    num_cpu = mp.cpu_count()
    pool = mp.Pool(num_cpu)
    pool.starmap_async(run_dock6, zip(pdb_ids))
    pool.close()
    pool.join()

    # pdb_id = '1AM0'
    # run_dock6(pdb_id)

    end = time.time()
    print('total elapsed time:', end - st, 'S')
