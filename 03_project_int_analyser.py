from My_Modules.structure.PDB import Protein
import os

current_dir = os.getcwd()

file_info = open("PDB_files/pdb_details")
src_dir = "PDB_files/src_pdb"
out_dir = "PDB_files/processed/"
for line in file_info:
    print(line)
    info = line.split()
    print(info)
    pdb_id = info[0] 
    print(pdb_id)
    pdb_name = pdb_id + ".pdb"
    print(pdb_name)
    ag_chain = info[1]
    ab_chain = info[2]
    

    antigen = Protein(os.path.join(src_dir,pdb_name),chain=ag_chain,format="PDB")
    antibody = Protein(os.path.join(src_dir,pdb_name),chain=ab_chain,format="PDB")
    complex = Protein(os.path.join(src_dir,pdb_name),chain=ag_chain+ab_chain,format="PDB")

    print(len(antigen.atoms),len(antibody.atoms),len(complex.atoms))

    dest_dir = os.path.join(out_dir,pdb_id)
    if not os.path.isdir(dest_dir):
        os.makedirs(dest_dir)
    antigen.write(dest_dir)
    antibody.write(dest_dir)
    complex.write(dest_dir)
    os.chdir(dest_dir)
    print(os.getcwd())
    antigen.run_naccess()
    antibody.run_naccess()
    complex.run_naccess()

    ag_asa = Protein(antigen.asa_file_name,format="ASA")
    ab_asa = Protein(antibody.asa_file_name,format="ASA")
    complex_asa = Protein(complex.asa_file_name,format="ASA")

    ag_asa.calc_asa_interface(complex_asa)
    ab_asa.calc_asa_interface(complex_asa)

    ag_buried_surface_area=0.0
    for atom in ag_asa.atoms:
        if atom.is_interface:
            ag_buried_surface_area += atom.asa_diff
    
    ab_buried_surface_area=0.0
    for atom in ab_asa.atoms:
        if atom.is_interface:
            ab_buried_surface_area += atom.asa_diff
    
    total_ag_asa = sum(atom.asa for atom in ag_asa.atoms)
    total_ab_asa = sum(atom.asa for atom in ab_asa.atoms)

    print("Total ASA of antigen:", total_ag_asa)
    print("Total ASA of antibody:", total_ab_asa)


            
    print(ag_buried_surface_area,ab_buried_surface_area)


    print(len(ag_asa.atoms),len(ab_asa.atoms),len(complex_asa.atoms))
    os.chdir(current_dir)


'''


    # antigen = PDB.parser(os.path.join(src_dir,pdb_name),chain = ag_chain)
    # antibody = PDB.parser(os.path.join(src_dir,pdb_name),chain = ab_chain)

    # complex = PDB.parser(os.path.join(src_dir,pdb_name))

     # print(len(antigen),len(antibody),len(complex))

    # ag_file_name = os.path.join(out_dir,pdb_id+"_"+ag_chain+".pdb")
    # PDB.write(antigen,ag_file_name)
    # ab_file_name = os.path.join(out_dir,pdb_id+"_"+ab_chain+".pdb")
    # PDB.write(antibody,ab_file_name)
    # complex_file_name = os.path.join(out_dir,pdb_id+"_"+ag_chain+ab_chain+".pdb") 
    # PDB.write(complex,complex_file_name)
    # os.chdir(out_dir)
    # os.system("naccess "+ag_file_name.split("/")[-1])
    # os.system("naccess "+ab_file_name.split("/")[-1])
    # os.system("naccess "+complex_file_name.split("/")[-1])

    ag_asa = PDB.parser(os.path.join(out_dir,pdb_id+"_"+ag_chain+".asa"))
    ab_asa = PDB.parser(os.path.join(out_dir,pdb_id+"_"+ab_chain+".asa"))
    complex_asa = PDB.parser(os.path.join(out_dir,pdb_id+"_"+ag_chain+ab_chain+".asa"))
    
    def check_n_update_asa(complex_atom,unbound_atom):
        if complex_atom == unbound_atom:
            if complex_atom.asa != unbound_atom.asa:
                unbound_atom.is_interface = True 
                unbound_atom.asa_cmplx = complex_atom.asa
                unbound_atom.asa_diff = unbound_atom.asa - complex_atom.asa
    
    for complex_atom in complex_asa:
        for ag_atom in ag_asa:
            check_n_update_asa(complex_atom,ag_atom)
            
        for ab_atom in ab_asa:
            check_n_update_asa(complex_atom,ab_atom)
    total_asa=0
    ag_buried_surface_area = 0
    outdir_2="PDB_results/outputs/03_project/naccess_result"
    if not os.path.isdir(outdir_2):
        os.makedirs(outdir_2)
    ag_int_output=open(os.path.join(outdir_2,"ag_naccess_result.pdb"),"w")
    for atom in ag_asa:
        if atom.is_interface:
            total_asa += atom.asa
            ag_buried_surface_area+=atom.asa_diff
            percentage = ag_buried_surface_area/total_asa
            out_str = "{0:}{1:>8.3f}{2:>8.3f}{3:>8.3f}\n".format(atom.source_line[:62],atom.asa_cmplx,atom.asa_diff,percentage)
            ag_int_output.write(out_str)

    # ag_asa = 
    # ab_asa = 
    # cmplx_asa = 

    # asa analysis code 

    # os.chdir(current_dir)
    '''    