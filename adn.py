import random


table = {
        'UUU': 'Phénylalanine', 'UUC': 'Phénylalanine',
        'UUA': 'Leucine', 'UUG': 'Leucine',
        'CUU': 'Leucine', 'CUC': 'Leucine', 'CUA': 'Leucine', 'CUG': 'Leucine',
        'AUU': 'Isoleucine', 'AUC': 'Isoleucine', 'AUA': 'Isoleucine',
        'AUG': 'Méthionine (Start)',
        'GUU': 'Valine', 'GUC': 'Valine', 'GUA': 'Valine', 'GUG': 'Valine',
        'UCU': 'Sérine', 'UCC': 'Sérine', 'UCA': 'Sérine', 'UCG': 'Sérine',
        'CCU': 'Proline', 'CCC': 'Proline', 'CCA': 'Proline', 'CCG': 'Proline',
        'ACU': 'Thréonine', 'ACC': 'Thréonine', 'ACA': 'Thréonine', 'ACG': 'Thréonine',
        'GCU': 'Alanine', 'GCC': 'Alanine', 'GCA': 'Alanine', 'GCG': 'Alanine',
        'UAU': 'Tyrosine', 'UAC': 'Tyrosine',
        'UAA': 'Stop', 'UAG': 'Stop',
        'CAU': 'Histidine', 'CAC': 'Histidine',
        'CAA': 'Glutamine', 'CAG': 'Glutamine',
        'AAU': 'Asparagine', 'AAC': 'Asparagine',
        'AAA': 'Lysine', 'AAG': 'Lysine',
        'GAU': 'Acide aspartique', 'GAC': 'Acide aspartique',
        'GAA': 'Acide glutamique', 'GAG': 'Acide glutamique',
        'UGU': 'Cystéine', 'UGC': 'Cystéine',
        'UGA': 'Stop',
        'UGG': 'Tryptophane',
        'CGU': 'Arginine', 'CGC': 'Arginine', 'CGA': 'Arginine', 'CGG': 'Arginine',
        'AGU': 'Sérine', 'AGC': 'Sérine',
        'AGA': 'Arginine', 'AGG': 'Arginine',
        'GGU': 'Glycine', 'GGC': 'Glycine', 'GGA': 'Glycine', 'GGG': 'Glycine'
    }



def generate_adn(length):
    bases = {
        'A': ['GCT', 'GCC', 'GCA', 'GCG'],
        'T': ['ACT', 'ACC', 'ACA', 'ACG'],
        'C': ['TGT', 'TGC', 'TGA', 'TGG'],
        'G': ['GGT', 'GGC', 'GGA', 'GGG', 'GAG', 'GAA', 'GAT', 'GAC']
    }
    random_adn_sequence = ''
    for _ in range(length):
        base = random.choice(list(bases.keys()))
        codon = random.choice(bases[base])
        random_adn_sequence += codon

    return random_adn_sequence


def read_adn_file(file_path):
    try:
        with open(file_path, 'r') as file:
            first_line = file.readline().strip()
            if first_line.startswith('>'):
                description = first_line[1:]
                adn_sequence = ''.join(line.strip() for line in file.readlines())
            else:
                adn_sequence = first_line
                for line in file:
                    adn_sequence += line.strip()

        return description, adn_sequence
    except FileNotFoundError:
        print(f"Le fichier '{file_path}' est introuvable.")
        return None, None


def is_adn_valid(adn_sequence):
    valid_bases = {'A', 'C', 'G', 'T'}
    if all(base in valid_bases for base in adn_sequence) :
        print("La séquence ADN est valide.")
    else:
        print("La séquence ADN est invalide.")    


def bases_freq(adn_sequence):
    base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    total_bases = len(adn_sequence)

    for base in adn_sequence:
        if base in base_counts:
            base_counts[base] += 1

    base_frequencies = {base: count / total_bases for base, count in base_counts.items()}
    return base_frequencies

def translate_to_arn(adn_sequence):
    return adn_sequence.replace('T', 'U')

def translate_to_proteins(arn_sequence):
    proteins = []
    for i in range(0, len(arn_sequence), 3):
        element = arn_sequence[i:i + 3]
        if element in table:
            amino_acid = table[element]
            if amino_acid == 'Stop':
                break 
            proteins.append(amino_acid)

    return proteins


def reverse_complement(adn_sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reversed_sequence = adn_sequence[::-1]
    complement_sequence = ''.join(complement[base] for base in reversed_sequence)
    return complement_sequence

def calculate_gc_content(adn_sequence):
    gc_count = adn_sequence.count('G') + adn_sequence.count('C')
    total_bases = len(adn_sequence)
    gc_content = (gc_count / total_bases) * 100
    return gc_content

def calculate_codon_frequencies(adn_sequence):
    codon_length = 3
    codon_counts = {}

    for i in range(0, len(adn_sequence) - codon_length + 1, codon_length):
        codon = adn_sequence[i:i + codon_length]
        if len(codon) == codon_length:
            if codon in codon_counts:
                codon_counts[codon] += 1
            else:
                codon_counts[codon] = 1

    total_codons = sum(codon_counts.values())
    codon_frequencies = {codon: count / total_codons for codon, count in codon_counts.items()}
    return codon_frequencies


def mutate_adn_sequence(adn_sequence, mutation_rate):
    mutation_types = ['A', 'C', 'G', 'T']
    mutated_sequence = ''

    for base in adn_sequence:
        if random.random() < mutation_rate:
            mutated_sequence += random.choice(mutation_types)
        else:
            mutated_sequence += base

    return mutated_sequence


def adn_motif(adn_sequence, motif):
    positions = []
    motif_length = len(motif)
    sequence_length = len(adn_sequence)

    for i in range(sequence_length - motif_length + 1):
        if adn_sequence[i:i + motif_length] == motif:
            positions.append(i)

    return positions

def generate_consensus_profile(aligned_sequences):
    profile_matrix = {'A': [], 'C': [], 'G': [], 'T': []}
    consensus_sequence = ''

    sequence_length = len(aligned_sequences[0])

    for i in range(sequence_length):
        column = [sequence[i] for sequence in aligned_sequences]
        for base in profile_matrix.keys():
            profile_matrix[base].append(column.count(base) / len(aligned_sequences))

        max_freq_base = max(profile_matrix.keys(), key=lambda x: profile_matrix[x][i])
        consensus_sequence += max_freq_base

    return consensus_sequence, profile_matrix
