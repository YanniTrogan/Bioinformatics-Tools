#DNA Toolset/Code testing file

from DNAToolkit import *
from Utilities import *
from structures import *
import random
from collections import *

print(f'DNAToolkit testing:\n')

randDNAStr = ''.join([random.choice(DNA_Nucleotides) for nuc in range(50)])

DNAStr = validateSeq(randDNAStr)

print(f"[4] + DNA String + Complement + Reverse Complement:\n5' {colored(DNAStr)} 3'")
print(f"   {''.join(['│' for c in range(len(DNAStr))])}")
print(f"3' {colored(complement(DNAStr))} 5'")
print(f"5' {colored(reverse_complement(DNAStr))} 3'")
print(f'[5] + GC Content:  {gc_content(DNAStr)}')
print(f'[6] + GC Content in Subsection k=5 {gc_content_subsec(DNAStr, k=5)}')
print(f'[7] + Aminoacids sequence from DNA: {translate_seq(DNAStr, 0)}\n')
print(f'[8] + Codon frequency (L): {codon_usage(DNAStr, "L")}')
    #finding the occurence of the specific codon "L"
print(f'[9] + Reading frames:')
for frame in gen_reading_frames(DNAStr):
    print(frame)

test_rf_frame = ['L', 'M', 'A', 'T', 'A', 'L', 'V', '_', 'V']
print(proteins_from_rf(test_rf_frame))

print(f'\n[10] + All proteins in 6 open reading frames:')
for prot in all_proteins_from_orfs(DNAStr, 0, 0, True):
    print(f'{prot}')

print(f'\nBio seq testing:\n')

from bio_seq import bio_seq
from bio_structs import *

test_dna = bio_seq("ATCG")

test_dna.get_seq_info()
test_dna.generate_rnd_seq(40, "RNA")

print(test_dna.get_seq_info())
print(test_dna.get_seq_biotype())
print(test_dna.nucleotide_frequency())
print(test_dna.transcription())
print(test_dna.reverse_complement())
print(test_dna.seq)
print(test_dna.gc_content())
print(test_dna.gc_content_subsec())
print(test_dna.translate_seq())
print(test_dna.codon_usage('L'))
for rf in test_dna.gen_reading_frames():
    print(rf)
print(test_dna.all_proteins_from_orfs())

writeTextFile("test.txt", test_dna.seq)
for rf in test_dna.gen_reading_frames():
    writeTextFile("test.txt", str(rf), 'a')