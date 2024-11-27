import argparse
from Bio import PDB
pdb_parser = PDB.PDBParser()

parser = argparse.ArgumentParser(
                    prog='All-Atom-Structure',
                    description='The program can change model from .pdb file in coarse-grained model to all-atom structure. Call a program with pdb file in coarse-grained model as a first argument.',
                    epilog='Program save result as a pdb file.')
parser.add_argument('pdb_file', help='input, pdb file in coarse-grained model')
parser.add_argument('--output', '-o', default = 'result.pdb', help='output, name of .pdb file')

args = parser.parse_args()
structure = pdb_parser.get_structure("cg_model", args.pdb_file)

newstructure = PDB.Structure.Structure("ResultStructure")
for model in structure:
  newmodel = PDB.Model.Model(model.get_id())
  newstructure.add(newmodel)
  for chain in model:
    newchain = PDB.Chain.Chain(chain.get_id())
    newmodel.add(newchain)
    for res in chain:
      newres = PDB.Residue.Residue(res.get_id(), res.get_resname(), res.get_segid())
      newchain.add(newres)

      atoms_rest = res.get_atoms()
      template_rest = []
      template_rest_less = []
      atoms_rest_less = []

      atoms = res.get_atoms()
      atoms_less = []
      template = []
      template_less = []

      match res.get_resname() :
        case "A" :
          struct_a = pdb_parser.get_structure("adenine", "struct_a2.pdb")
          for a in struct_a.get_atoms():
            template.append(a)
          for a in template:
            if a.get_name() in ['N9', 'C2', 'C6']:
              template_less.append(a);
          for a in atoms:
            if a.get_name() in ['N9', 'C2', 'C6']:
              atoms_less.append(a)

          struct_rest = pdb_parser.get_structure("rest", "struct_rest_purines.pdb")
          for a in struct_rest.get_atoms():
            template_rest.append(a)
          for a in template_rest:
            if a.get_name() in ['P', 'C4\'', 'N9']:
              template_rest_less.append(a);
          for a in atoms_rest:
            if a.get_name() in ['P', 'C4\'', 'N9']:
              atoms_rest_less.append(a)

        case "G" :
          struct_g = pdb_parser.get_structure("guanine", "struct_g2.pdb")
          for a in struct_g.get_atoms():
            template.append(a)
          for a in template:
            if a.get_name() in ['N9', 'C2', 'C6']:
              template_less.append(a);
          for a in atoms:
            if a.get_name() in ['N9', 'C2', 'C6']:
              atoms_less.append(a)

          struct_rest = pdb_parser.get_structure("rest", "struct_rest_purines.pdb")
          for a in struct_rest.get_atoms():
            template_rest.append(a)
          for a in template_rest:
                if a.get_name() in ['P', 'C4\'', 'N9']:
                  template_rest_less.append(a);
          for a in atoms_rest:
                if a.get_name() in ['P', 'C4\'', 'N9']:
                  atoms_rest_less.append(a)

        case "C" :
          struct_c = pdb_parser.get_structure("cytosine", "struct_c2.pdb")
          for a in struct_c.get_atoms():
            template.append(a)
          for a in template:
            if a.get_name() in ['N1', 'C2', 'C4']:
              template_less.append(a);
          for a in atoms:
            if a.get_name() in ['N1', 'C2', 'C4']:
              atoms_less.append(a)

          struct_rest = pdb_parser.get_structure("rest", "struct_rest_pyrimidines.pdb")
          for a in struct_rest.get_atoms():
            template_rest.append(a)
          for a in template_rest:
                if a.get_name() in ['P', 'C4\'', 'N1']:
                  template_rest_less.append(a);
          for a in atoms_rest:
                if a.get_name() in ['P', 'C4\'', 'N1']:
                  atoms_rest_less.append(a)

        case "U" :
          struct_u = pdb_parser.get_structure("uracil", "struct_u2.pdb")
          for a in struct_u.get_atoms():
            template.append(a)
          for a in template:
            if a.get_name() in ['N1', 'C2', 'C4']:
              template_less.append(a);
          for a in atoms:
            if a.get_name() in ['N1', 'C2', 'C4']:
              atoms_less.append(a)

          struct_rest = pdb_parser.get_structure("rest", "struct_rest_pyrimidines.pdb")
          for a in struct_rest.get_atoms():
            template_rest.append(a)
          for a in template_rest:
                if a.get_name() in ['P', 'C4\'', 'N1']:
                  template_rest_less.append(a);
          for a in atoms_rest:
                if a.get_name() in ['P', 'C4\'', 'N1']:
                  atoms_rest_less.append(a)
        case _ :
          break
      
      sup = PDB.Superimposer()
      sup.set_atoms(fixed = atoms_less, moving=template_less)
      sup.apply(template)
      for a in template:
        newres.add(a)

      everything_is_ok = 0
      for a in atoms_rest_less:
        if a.get_name() == 'P':
          everything_is_ok = 1
      if everything_is_ok == 0:
        template_rest_less_copy = template_rest_less.copy()
        template_rest_less.clear()
        for a in template_rest_less_copy:
          if a.get_name() != 'P':
            template_rest_less.append(a)

      sup_rest = PDB.Superimposer()
      sup_rest.set_atoms(fixed = atoms_rest_less, moving=template_rest_less)
      template_rest2 = []
      for a in template_rest:
        if a.get_name() != 'N9' and a.get_name() != 'N1':
          template_rest2.append(a)
      sup_rest.apply(template_rest2)
      for a in template_rest2:
        newres.add(a)


io = PDB.PDBIO()
io.set_structure(newstructure)
io.save(args.output)
print("\nThe result has been saved to the "+args.output+" file.")
