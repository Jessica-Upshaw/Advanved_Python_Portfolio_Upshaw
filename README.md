# Advanved_Python_Portfolio_Upshaw
This is the portfolio of python code that I learned during BISC 450C 

```python
# Part 1
```


```python
# Identifying a sequence to work with
```


```python
# pip install biopython in the terminal
from Bio.Seq import Seq
```


```python
# Run code with shift+enter
```


```python
# my_seq is a variable that represents a sequence
my_seq = Seq("GATCG")
```


```python
for index, letter in enumerate(my_seq):
    print("%i %s" % (index, letter))
```

    0 G
    1 A
    2 T
    3 C
    4 G



```python
print(len(my_seq))
```

    5



```python
print(my_seq[0])
```

    G



```python
print(my_seq[4])
```

    G



```python
# Counting how many time AA occurs
# Doesn't count overlapping
Seq("AAAA").count("AA")
```




    2




```python
my_seq = Seq("GAGACGTCAGCTACGACTACGATCAGCAT")
```


```python
len(my_seq)
```




    29




```python
my_seq.count("G")
```




    7




```python
# Calculating G+C content 
100 * (my_seq.count("G") + my_seq.count("C")) / len(my_seq)
```




    51.724137931034484




```python
from Bio.SeqUtils import gc_fraction
```


```python
my_seq = Seq("GAGACGTCAGCTACGACTACGATCAGCAT")
```


```python
# Easier way to caculate G+C content
gc_fraction(my_seq)
```




    0.5172413793103449




```python
# Splicing nucleotides
```


```python
# [start:stop:step]
# Start at 0, get every 3rd nucleotide
# Start with G, skip AG, A, skip CG, T, etc
my_seq[0::3]
```




    Seq('GATGAAAAAA')




```python
# Empty start, starts at the beginning by default
my_seq[::3]
```




    Seq('GATGAAAAAA')




```python
my_seq[1::3]
```




    Seq('ACCCCCCTGT')




```python
my_seq[2::6]
```




    Seq('GAGGC')




```python
# [start:stop]
# Start at the start element
# Will include elements up to but not including the stop element
# All elements starting at 2 and up to but not including 3
my_seq[2:3]
```




    Seq('G')




```python
# Use negative to start from the other end
# Just prints the sequence backwards 
my_seq[::-1]
```




    Seq('TACGACTAGCATCAGCATCGACTGCAGAG')




```python
# Start at the end of the sequence, get every 3rd element, working backwards
my_seq[::-3]
```




    Seq('TGTCCCCCCA')




```python
# Turning seq object back into string
str(my_seq)
```




    'GAGACGTCAGCTACGACTACGATCAGCAT'




```python
# Gave it a name
# \n is a newline character that moves the cursor to the next line
# %s is a placeholder for the string identified with %
fasta_format_string = ">Tom\n%s\n" % my_seq
```


```python
print(fasta_format_string)
```

    >Tom
    GAGACGTCAGCTACGACTACGATCAGCAT
    



```python
# Took out the \n
# Everything on the same line
fasta_format_string = ">Tom%s\n" % my_seq
```


```python
print(fasta_format_string)
```

    >TomGAGACGTCAGCTACGACTACGATCAGCAT
    



```python
# Didn't identify what %s is a placeholder for
fasta_format_string = ">Tom\n%s\n"
```


```python
print(fasta_format_string)
```

    >Tom
    %s
    



```python
# Adding two scripts together
seq1 = Seq("AGTC")
seq2 = Seq("GTCA")
```


```python
seq1 + seq2
```




    Seq('AGTCGTCA')




```python
seq2 + seq1
```




    Seq('GTCAAGTC')




```python
# Part 2
```


```python
contigs = [Seq("ATG"), Seq("ATCCCG"), Seq("TTTGCA")]
```


```python
# Sequence of ten unknown nucleotides
spacer = Seq("N" *10)
```


```python
# Take the spacer object that we made and join with the contigs
# ATG-N*10-ATCCCG-N*10-TTTGCA
spacer.join(contigs)
```




    Seq('ATGNNNNNNNNNNATCCCGNNNNNNNNNNTTTGCA')




```python
# Case sensitivity
dna_seq = Seq("acgtACGT")
```


```python
dna_seq
```




    Seq('acgtACGT')




```python
dna_seq.upper()
```




    Seq('ACGTACGT')




```python
dna_seq.lower()
```




    Seq('acgtacgt')




```python
# Save the sequence as uppercase
dna_seq = dna_seq.upper()
```


```python
# False because it's case sensitive
# The sequence has gtAC not gtac
"gtac" in dna_seq
```




    False




```python
"GTAC" in dna_seq
```




    True




```python
my_seq = Seq("GATCGATCGATCTCGATCGATCGATCGATCGATCGTT")
```


```python
my_seq.complement()
```




    Seq('CTAGCTAGCTAGAGCTAGCTAGCTAGCTAGCTAGCAA')




```python
my_seq.reverse_complement()
```




    Seq('AACGATCGATCGATCGATCGATCGAGATCGATCGATC')




```python
# Ambiguity codes
# R: A, G
# Y: C, T
# S: G, C
# W: A, T
# K: G, T
# M: A, C
# B: C, G, T
# D: A, G, T
# H: A, C, T
# V: A, C, G
# N: A, C, G, T
protein_seq = Seq("EVRNAK")
protein_seq.complement()
```




    Seq('EBYNTM')




```python
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
```


```python
coding_dna
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
# Template strand is the original dna
# Coding strand is the reverse-complement of template
# mRNA is the transcribed coding strand
template_dna = coding_dna.reverse_complement()
```


```python
template_dna
```




    Seq('CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT')




```python
coding_dna
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
messenger_rna = coding_dna.transcribe()
```


```python
messenger_rna
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
# Get the same thing
template_dna.reverse_complement().transcribe()
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
# Reverses transcription
# Back to coding dna
messenger_rna.back_transcribe()
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
messenger_rna
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
# Astericks are stop codons
messenger_rna.translate()
```




    Seq('MAIVMGR*KGAR*')




```python
# Part 3
```


```python
# Using mitochondrial codon table
coding_dna.translate(table="Vertebrate Mitochondrial")
```




    Seq('MAIVMGRWKGAR*')




```python
coding_dna.translate(table = 2)
```




    Seq('MAIVMGRWKGAR*')




```python
# Stops at the stop codon
coding_dna.translate(to_stop = True)
```




    Seq('MAIVMGR')




```python
coding_dna.translate(table = 2, to_stop = True)
```




    Seq('MAIVMGRWKGAR')




```python
coding_dna.translate(table = 2, stop_symbol = "!")
```




    Seq('MAIVMGRWKGAR!')




```python
coding_dna.translate(table = 2, stop_symbol = "<3")
```




    Seq('MAIVMGRWKGAR<3')




```python
# Non standard stop codon
gene = Seq("GTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCAGCACAGGCTGCGGAAATTACGTTAGTCCCGTCAGTAAAATTACAGATAGGCGATCGTGATAATCGTGGCTATTACTGGGATGGAGGTCACTGGCGCGACCACGGCTGGTGGAAACAACATTATGAATGGCGAGGCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCATAAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAA")
```


```python
gene.translate(table = "Bacterial")
```




    Seq('VKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HR*')




```python
gene.translate(table = "Bacterial", to_stop = True)
```




    Seq('VKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HHR')




```python
gene.translate(table = "Bacterial", cds = True)
```




    Seq('MKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HHR')




```python
from Bio.Data import CodonTable
```


```python
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
```


```python
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
```


```python
print(standard_table)
```

    Table 1 Standard, SGC0
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
    T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I   | ACT T   | AAT N   | AGT S   | T
    A | ATC I   | ACC T   | AAC N   | AGC S   | C
    A | ATA I   | ACA T   | AAA K   | AGA R   | A
    A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V   | GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
print(mito_table)
```

    Table 2 Vertebrate Mitochondrial, SGC1
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA W   | A
    T | TTG L   | TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L   | CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I(s)| ACT T   | AAT N   | AGT S   | T
    A | ATC I(s)| ACC T   | AAC N   | AGC S   | C
    A | ATA M(s)| ACA T   | AAA K   | AGA Stop| A
    A | ATG M(s)| ACG T   | AAG K   | AGG Stop| G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V(s)| GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
mito_table.stop_codons
```




    ['TAA', 'TAG', 'AGA', 'AGG']




```python
mito_table.start_codons
```




    ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']




```python
# Sequence comparison
seq = Seq("AGTC")
```


```python
# These are fully equal
"AGTC" == seq1
```




    True




```python
seq1 == "AGTC"
```




    True




```python
# Length of sequence known but not the letters
unknown_seq = Seq(None, 10)
```


```python
unknown_seq
```




    Seq(None, length=10)




```python
len(unknown_seq)
```




    10




```python
# Part 4
```


```python
seq = Seq({117512683: "TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT"}, length = 159345973)
```


```python
# Undefined
seq[1000:1020]
```




    Seq(None, length=20)




```python
# Defined
seq[117512690:117512700]
```




    Seq('CCTGAATGTG')




```python
seq[117512670:]
```




    Seq({13: 'TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT'}, length=41833303)




```python
seq = Seq("AGTC")
```


```python
undefined_seq = Seq(None, length = 10)
```


```python
seq + undefined_seq + seq
```




    Seq({0: 'AGTC', 14: 'AGTC'}, length=18)




```python
my_seq = Seq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
```


```python
from Bio.Seq import MutableSeq
```


```python
mutable_seq = MutableSeq(my_seq)
```


```python
mutable_seq
```




    MutableSeq('GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
# Fifth position (technically sixth) is now C
mutable_seq[5] = "C"
```


```python
mutable_seq
```




    MutableSeq('GCCATCGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
# Only removes the first T
mutable_seq.remove("T")
```


```python
mutable_seq
```




    MutableSeq('GCCACGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
mutable_seq.reverse()
```


```python
mutable_seq
```




    MutableSeq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')




```python
# If you want your new code to be free from mutation/change
new_seq = Seq(mutable_seq)
```


```python
# Back to a normal seq
new_seq
```




    Seq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')




```python
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
```


```python
my_string = "GTAGCTAGCTAGATCGATCGATCGAGGTACAGGATATCG"
```


```python
reverse_complement(my_string)
```




    'CGATATCCTGTACCTCGATCGATCGATCTAGCTAGCTAC'




```python
transcribe(my_string)
```




    'GUAGCUAGCUAGAUCGAUCGAUCGAGGUACAGGAUAUCG'




```python
translate(my_string)
```




    'VAS*IDRSRYRIS'




```python
back_transcribe(my_string)
```




    'GTAGCTAGCTAGATCGATCGATCGAGGTACAGGATATCG'




```python
# Done
```
