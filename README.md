# Coarse-grained-model
Dwa skrypty pisane na potrzeby przedmiotu bioinformatyka strukturalna. Pierwszy przekształca strukturę RNA zawartą w pliku .pdb na model gruboziarnisty, a drugi odwrotnie, próbuje wrócić z modelu gruboziarnistego do pełnoatomowej struktury DNA.

Żeby drugi skrypt (do struktury pełnoatomowej) działał, potrzebne są ściągnięte szablony struktur (pliki pdb zaczynające się od struct_) umiejscowione w tym samym katalogu co skrypt.

# coarse_grained.py
Przykładowe uruchomienie skryptu z lini poleceń:

python coarse_grained.py 430D --output_dir ./pdb_files --output_file 430D_coarse_grained.pdb

Argumenty:

   pdb_code: wymagany, kod struktury PDB (np. 430D).
   --output_dir: opcjonalny, katalog docelowy na plik PDB i wynik (domyślnie bieżący katalog).
   --output_file: opcjonalny, nazwa pliku wyjściowego (domyślnie coarse_grained.pdb)
