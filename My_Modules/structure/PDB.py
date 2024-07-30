import os
class Atom():
    def __init__(self,atom_line):
        self.x_coord = float(atom_line[30:38])
        self.y_coord = float(atom_line[38:46])
        self.z_coord = float(atom_line[46:54])
        
        self.chain = atom_line[21]
        
        self.name=atom_line[12:16].strip()
        self.number=atom_line[6:11].strip()
        
        self.residue_name=atom_line[17:20].strip()
        self.residue_number=atom_line[22:26].strip()
        
        self.source_line=atom_line
            
        self.is_interface=False

    
    def __eq__(self,other):
        if isinstance(other,Atom):
            return (self.name == other.name and
                    self.number == other.number and
                     self.residue_name == other.residue_name and
                      self.residue_number == other.residue_number and
                       self.chain == other.chain )
        

class Atom_asa(Atom):
    def __init__(self, atom_line):
        super().__init__(atom_line)
        self.asa = float(atom_line[54:62])
        self.asa_complex = None
        self.asa_diff = None
    
                
class Protein():
    def __init__(self,file_path,chain=None,format=None):
        self.file_path = file_path
        self.file_name = self.file_path.split("/")[-1]
        self.id = self.file_name.split(".")[0]
        self.extension = self.file_name.split(".")[-1]
        self.chain = chain 
        if self.chain:
            self.custom_name = self.id+"_"+self.chain+"."+self.extension
            self.custom_id = self.id+"_"+self.chain
        else:
            self.custom_name = self.id+"."+self.extension
            self.custom_id = self.id
        self.format = format 
        self.atoms = []
        self.residues = []

        self.parse()

    def parse(self):
        print("File path is: ",self.file_path)

        source = open(self.file_path).readlines()
        print("Chain provided is: ",self.chain)
        
        for line in source:
            if line.startswith("ATOM"):
                if self.chain and line[21] in self.chain:
                    if self.format == "PDB":
                        self.atoms.append(Atom(line))
                    if self.format == "ASA":
                        self.atoms.append(Atom_asa(line))
                    
                
                if self.chain is None:   # if not chain
                    if self.format == "PDB":
                        self.atoms.append(Atom(line))
                    if self.format == "ASA":
                        self.atoms.append(Atom_asa(line))

    def run_naccess(self):
        self.asa_file_name = self.custom_id + ".asa"
        os.system("naccess "+self.custom_name)

    def write(self,dest_dir):
        self.out_file_name = os.path.join(dest_dir,self.custom_name)
        outfile = open(self.out_file_name,'w')
        for atom in self.atoms:
            outfile.write(atom.source_line)
                

    def calc_asa_interface(self,other):
        if isinstance(other,Protein):
            for atom_a in self.atoms:
                for complex_atom in other.atoms:
                    if atom_a == complex_atom and atom_a.asa!=complex_atom.asa:
                        atom_a.is_interface = True 
                        complex_atom.is_interface = True 
                        atom_a.asa_diff = atom_a.asa - complex_atom.asa 
                        atom_a.asa_complex = complex_atom.asa 


