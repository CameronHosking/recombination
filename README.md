# recombination
detects recombination events
# usage
python3 find-mutations-fast.py "210118trimmed.noDup.mafftAlignment" > "210118mutations.txt"
./same-mutations "210118mutations.txt" "210118.recombinations2" "210118.duplicates.txt"
