with open('rosalind_gc.txt') as f:
    fasta_data = f.readlines()
def parse_fasta(fasta_data):
    sequences = {}
    current_id = None
    for line in fasta_data:
        line = line.strip()
        if line.startswith('>'):
            current_id = line[1:]
            sequences[current_id] = ''
        else:
            sequences[current_id] += line
    return sequences
sequences = parse_fasta(fasta_data)
max_gc_content = 0
max_gc_id = None
for seq_id, sequence in sequences.items():
    gc_count = sequence.count('G') + sequence.count('C')
    gc_content = (gc_count / len(sequence)) * 100
    if gc_content > max_gc_content:
        max_gc_content = gc_content
        max_gc_id = seq_id
print(max_gc_id)
print(max_gc_content)