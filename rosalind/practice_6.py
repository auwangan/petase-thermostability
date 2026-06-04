with open ('rosalind_hamm.txt','r') as file:
    lines = file.readlines()
s = lines[0].strip()
t = lines[1].strip()
def calculate_hamming_distance(s, t):
    count =0
    for i in range(len(s)):
        if s[i] != t[i]:
            count += 1 
    return count
print(calculate_hamming_distance(s, t))
def calculate_hamming_distance(s, t):