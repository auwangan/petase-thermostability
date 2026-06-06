with open("rosalind_revc.txt") as f:
    DNA_SEQUENCE = f.read().strip()
DNA_COMPLEMENT = DNA_SEQUENCE.replace("A","t").replace("T","a").replace("C","g").replace("G","c").upper()
DNA_REVERSE_COMPLEMENT = DNA_COMPLEMENT[::-1]
print(DNA_REVERSE_COMPLEMENT)