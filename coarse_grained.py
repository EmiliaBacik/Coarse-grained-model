import argparse
from Bio import PDB
from Bio.PDB import PDBList, PDBParser, PDBIO, Select

# Parser argumentów wiersza poleceń
parser = argparse.ArgumentParser(description="Opis działania skryptu.")
parser.add_argument("pdb_code", help="Kod PDB struktury do pobrania (np. '430D').")
parser.add_argument("--output_dir", default=".", help="Katalog docelowy dla pobranego pliku PDB.")
parser.add_argument("--output_file", default="coarse_grained.pdb", help="Nazwa pliku wyjściowego.")
args = parser.parse_args()

# Pobieranie struktury PDB
pdbl = PDBList()
pdb_file = pdbl.retrieve_pdb_file(args.pdb_code, file_format='pdb', pdir=args.output_dir)
print(f"Pobrany plik PDB: {pdb_file}")

pdb_parser = PDBParser(QUIET=True)
structure = pdb_parser.get_structure(args.pdb_code, pdb_file)

# Klasa do filtracji atomów
class CoarseGrainedSelect(Select):

    def accept_atom(self, atom):
        purines_atoms = {"N9", "C2", "C6"}
        pyrimidines_atoms = {"N1", "C2", "C4"}
        backbone_atoms = {"P", "C4'"}

        resname = atom.get_parent().get_resname()
        atom_name = atom.get_name()

        if resname in {"A", "G"}:
            return atom_name in purines_atoms or atom_name in backbone_atoms
        elif resname in {"C", "U"}:
            return atom_name in pyrimidines_atoms or atom_name in backbone_atoms
        return False

# Zapisanie reprezentacji gruboziarnistej do pliku
output_path = f"{args.output_dir}/{args.output_file}"
io = PDBIO()
io.set_structure(structure)
io.save(output_path, select=CoarseGrainedSelect())
print(f"Reprezentacja gruboziarnista zapisana do pliku: {output_path}")



