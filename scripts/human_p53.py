import random
import os

RNA_codon_table = {}

script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.normpath(os.path.join(script_dir, '..', 'data'))
results_dir = os.path.normpath(os.path.join(script_dir, '..', 'results'))

#### functions ####

# function that replaces a nucleotide
def Mutate_DNA(seq):
    if len(seq) == 0:
        return seq
    nuc_list = ["A", "T", "G", "C"]
    seq = seq.upper()

    rand_index = random.randrange(len(seq))
    rand_nuc = seq[rand_index]
    nuc_list.remove(rand_nuc)
    new_nuc = random.choice(nuc_list)

    return seq[:rand_index] + new_nuc + seq[rand_index + 1:]

# function that adds a nucleotide
def Insert_DNA(seq):
    if len(seq) == 0:
        return seq
    nuc_list = ["A", "T", "G", "C"]
    seq = seq.upper()

    rand_index = random.randrange(len(seq) + 1)
    new_nuc1 = random.choice(nuc_list)
    new_nuc2 = random.choice(nuc_list)
    new_nucs = new_nuc1 + new_nuc2
    return seq[:rand_index] + new_nucs + seq[rand_index:]

# function that deletes a nucleotide
def Delete_DNA(seq):
    if len(seq) == 0:
        return seq
    seq = seq.upper()

    k = random.randint(1, 3)
    rand_index = random.randrange(len(seq))
    if k == 1:
        return seq[:rand_index] + seq[rand_index + 1:]
    else:
        return seq[:rand_index] + seq[rand_index + k:]

# compare sequences
def Comp_seq(a, b):
    count = abs(len(a) - len(b))
    for x, y in zip(a, b):
        if x != y:
            count += 1
    return count

# reads codon table
def Read_dict():
    global RNA_codon_table
    with open(os.path.join(data_dir, "codon_AA.txt")) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                codon = parts[0].strip()
                aa = parts[1].strip()
                RNA_codon_table[codon] = aa

# transcribes DNA to RNA
def DNA_RNA_Cod(dna_seq):
    return dna_seq.upper().replace("T", "U")

# translates RNA to protein
def RNA_prot(rna_seq):
    protein = ""
    for i in range(0, len(rna_seq) - 2, 3):
        codon = rna_seq[i:i+3]
        if codon not in RNA_codon_table:
            break
        aa = RNA_codon_table[codon]
        protein += aa
    return protein

#### main code #####

def main():
    Read_dict()

    # read p53 DNA
    orgnl_seq = ""
    with open(os.path.join(data_dir, "human_p53_coding.txt"), "r") as f:
        for line in f:
            line = line.strip().upper()
            if line.startswith(">"):
                continue
            orgnl_seq += line

    orig_prot = RNA_prot(DNA_RNA_Cod(orgnl_seq))

    runs = 1000
    total_generations = 0

    for run in range(runs):
        generations = 0
        mutated_seq = orgnl_seq

        while True:
            generations += 1

            # מוטציה מתרחשת רק בהסתברות של 10^-4
            if random.random() < 0.0001:

                mutations = random.random()

                # בחירת סוג מוטציה
                if mutations < 0.98:
                    mutated_seq = Mutate_DNA(mutated_seq)
                elif mutations < 0.99:
                    mutated_seq = Insert_DNA(mutated_seq)
                    break
                else:
                    mutated_seq = Delete_DNA(mutated_seq)
                    break

                # תרגום לחלבון אחרי המוטציה
                mutated_prot = RNA_prot(DNA_RNA_Cod(mutated_seq))

                # אם החלבון השתנה – עוצרים את הלולאה
                if Comp_seq(orig_prot, mutated_prot) > 0:
                    break

        total_generations += generations

    avg_generations = total_generations / runs

    with open(os.path.join(results_dir, "results.txt"), "w") as f:
        f.write("Average generations until mutation changes P53 protein: " + str(avg_generations) + "\n")


main()
