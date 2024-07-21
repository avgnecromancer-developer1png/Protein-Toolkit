
List=['A', 'T', 'G', 'C', 'R', 'N', 'D', 'E', 'Q', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'W', 'Y', 'V']

# Validate Protein Sequence
def validateSeq(aa_seq):
   AAseq=aa_seq.upper()
   for aa in AAseq:
        if aa not in List:
             return False
   return AAseq


randseq = input("Enter sequence:\n")
print("\nProtein sequence is:\n",validateSeq(randseq))

#Lengh of Protein
print("\nlength of protein sequence is:")
print(len (validateSeq(randseq)))

# Frequency of Aminoacids in Protein
def countAmiFrequency(seq):
    seq_upper=seq.upper()
    SeqFreq = {'A':0, 'T':0, 'G':0, 'C':0, 'R':0, 'N':0, 'D':0, 'E':0, 'Q':0, 'H':0, 'I':0, 'L':0, 'K':0, 'M':0, 'F':0, 'P':0, 'S':0, 'W':0, 'Y':0, 'V':0}
    for n in seq_upper:
        SeqFreq[n] += 1
    return SeqFreq
print("\nFrequency of each Amino Acid is\n",countAmiFrequency(randseq))

# Protein -> RNA
def RNA_seq(seq):
    seq_upper=seq.upper()
    RNA_Codons = {
        # 'M' - START, '_' - STOP
        "A": "GCU", "A": "GCC", "A": "GCA", "A": "GCG",
        "C": "UGU", "C": "UGC",
        "D": "GAU", "D": "GAC",
        "E": "GAA", "E": "GAG",
        "F": "UUU", "F": "UUC",
        "G": "GGU", "G": "GGC", "G": "GGA", "G": "GGG",
        "H": "CAU", "H": "CAC",
        "I": "AUA", "I": "AUU", "I": "AUC",
        "K": "AAA", "K": "AAG",
        "L": "UUA", "L": "UUG", "L": "CUU", "L": "CUC", "L": "CUA", "L": "CUG",
        "M": "AUG",
        "N": "AAU", "N": "AAC",
        "P": "CCU", "P": "CCC", "P": "CCA", "P": "CCG",
        "Q": "CAA", "Q": "CAG",
        "R": "CGU", "R": "CGC", "R": "CGA", "R": "CGG", "R": "AGA", "R": "AGG",
        "S": "UCU", "S": "UCC", "UCA": "S", "TCG": "S", "AGU": "S", "AGC": "S",
        "T": "ACU", "T": "ACC", "T": "ACA", "T": "ACG",
        "V": "GUU", "V": "GUC", "GUA": "V", "GUG": "V",
        "W": "UGG",
        "Y": "UAU", "Y": "UAC",
        "_": "UAA", "_": "UAG", "_": "TGA"
    }

    return ''.join([RNA_Codons[seq_upper[pos:pos + 1]] for pos in range( len(seq_upper))])


print("\nRNA Sequence is:\n",RNA_seq(randseq))

# RNA -> DNA
def DNA_seq(seq):
    seq_upper=seq.upper()
    DNA_Codons = {
        # 'M' - START, '_' - STOP
        "A": "GCT", "A": "GCC", "A": "GCA", "A": "GCG",
        "C": "TGT", "C": "TGC",
        "D": "GAT", "D": "GAC",
        "E": "GAA", "E": "GAG",
        "F": "TTT", "F": "TTC",
        "G": "GGT", "G": "GGC", "G": "GGA", "G": "GGG",
        "H": "CAT", "H": "CAC",
        "I": "ATA", "I": "ATT", "I": "ATC",
        "K": "AAA", "K": "AAG",
        "L": "TTA", "L": "TTG", "L": "CTT", "L": "CTC", "L": "CTA", "L": "CTG",
        "M": "ATG",
        "N": "AAT", "N": "AAC",
        "P": "CCT", "P": "CCC", "P": "CCA", "P": "CCG",
        "Q": "CAA", "Q": "CAG",
        "R": "CGT", "R": "CGC", "R": "CGA", "R": "CGG", "R": "AGA", "R": "AGG",
        "S": "TCT", "S": "TCC", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
        "T": "ACT", "T": "ACC", "T": "ACA", "T": "ACG",
        "V": "GTT", "V": "GTC", "GTA": "V", "GTG": "V",
        "W": "TGG",
        "Y": "TAT", "Y": "TAC",
        "_": "TAA", "_": "TAG", "_": "TGA"
    }
    return ''.join([DNA_Codons[seq_upper[pos:pos + 1]] for pos in range(len(seq_upper))])

#(DNA_seq(randseq))
print("\nDNA sequence is:\n",(DNA_seq(randseq)))

# Complement
print("\nComplement of DNA sequence is:")
Comp_DNAseq =DNA_seq(randseq)

for a in Comp_DNAseq:
    if (a == "A"):
        print("T",end='')
    elif (a == "T"):
            print("A",end='')
    elif (a == "G"):
            print("C",end='')
    elif (a == "C"):
            print("G",end='')

# Reverse Complement
print("\n\nReverse Complement of DNA sequence is:")
Reverse_CompDNASeq=Comp_DNAseq[::-1]
for a in Reverse_CompDNASeq:
    if (a == "A"):
        print("T",end='')
    elif (a == "T"):
            print("A",end='')
    elif (a == "G"):
            print("C",end='')
    elif (a == "C"):
            print("G",end='')

# GC content
def gc_Content(seq):
  ASeq=seq.upper()
  return round((ASeq.count("G") + ASeq.count("C"))/len(ASeq)* 100)
print("\n\nGC% content is:\n",gc_Content(DNA_seq(randseq)),"%\n")

input("Please press enter to close.")



