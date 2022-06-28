from bio_seq import bio_seq
from utilities import readFASTA, readTextFile, writeTextFile

test_DNA = bio_seq(sequence, label= label_name)
print(test_DNA.get_seq_info())

#test_DNA.generate_rnd_seq(40, "RNA")
#print(test_DNA.get_seq_info())

print(test_DNA.nuc_freq())
print(test_DNA.transcription())
print(test_DNA.seq)
print(test_DNA.reverse_complement())
print(test_DNA.gc_content())
print(test_DNA.gc_content_subseq())
print(test_DNA.translate_seq())
print(test_DNA.codon_freq("L"))
for rf in test_DNA.gen_reading_frames():
    print(rf)
#print(test_DNA.proteins_from_rf(['G', 'G', 'M', 'K', 'A', '_', 'L', 'D', 'F', 'L', 'A', '_', 'L']))
print(test_DNA.all_prots_from_orfs())

#writeTextFile("test.txt", test_DNA.seq)
#for rf in test_DNA.gen_reading_frames():
#    writeTextFile("test.txt", str(rf), 'a')






