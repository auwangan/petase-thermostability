with open("rosalind_dna.txt") as f:
    DNA_SEQUENCE = f.read().strip()
COUNT_A = DNA_SEQUENCE.count("A")
COUNT_C = DNA_SEQUENCE.count("C")
COUNT_G = DNA_SEQUENCE.count("G")
COUNT_T = DNA_SEQUENCE.count("T")   
print(COUNT_A, COUNT_C, COUNT_G, COUNT_T)