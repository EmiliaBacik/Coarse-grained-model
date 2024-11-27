# Coarse-grained-model
Dwa skrypty pisane na potrzeby przedmiotu bioinformatyka strukturalna. Pierwszy przekształca strukturę RNA zawartą w pliku .pdb na model gruboziarnisty, a drugi odwrotnie, próbuje wrócić z modelu gruboziarnistego do pełnoatomowej struktury DNA.

UWAGA!
Do działania skryptów potrzebna jest zainstalowana biblioteka biopython. 
Żeby drugi skrypt (do struktury pełnoatomowej) działał, potrzebne są ściągnięte szablony struktur (pliki pdb zaczynające się od struct_) umiejscowione w tym samym katalogu co skrypt.

# coarse_grained.py
Przekształcenie pełnoatomowej struktury RNA na model gruboziarnisty.
Przykładowe uruchomienie skryptu z lini poleceń:

python coarse_grained.py 430D --output_dir ./pdb_files --output_file 430D_coarse_grained.pdb

Argumenty:

   pdb_code: wymagany, kod struktury PDB (np. 430D).
   --output_dir: opcjonalny, katalog docelowy na plik PDB i wynik (domyślnie bieżący katalog).
   --output_file: opcjonalny, nazwa pliku wyjściowego (domyślnie coarse_grained.pdb)


# all-atom.py
Przekształcenie modelu gruboziarnistego RNA na pełnoatomową strukturę.
Przykładowe uruchomienie skryptu z lini poleceń:

python all_atom.py 430D_coarse_grained.pdb --output_file 430D_full_atom.pdb

Argumenty:

   pdb_file: wymagany, plik z modelem gruboziarnistym struktury RNA zapisany w formacie PDB.
   --output_file: opcjonalny, nazwa pliku wyjściowego (domyślnie result.pdb)
