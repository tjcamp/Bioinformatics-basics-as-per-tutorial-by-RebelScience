import random
import collections
from collections import Counter
from bio_structs import DNA_Codons, NUCLEOTIDE_BASE, RNA_Codons


class bio_seq:
    """DNA sequence class. Default value: ATCG, DNA, No label"""

    def __init__(self, seq="ATCG", seq_type="DNA", label="No Label"):
        """Sequence initialization validation"""
        self.seq = seq.upper()
        self.label = label
        self.seq_type = seq_type
        self.is_valid = self.__validate()
        assert self.is_valid, f"Provided data does not seem to be a correct {self.seq_type} sequence"

    # DNA toolkit functions
    def __validate(self):
        """To show if the input sequence is a DNA seq or not"""
        return set(NUCLEOTIDE_BASE[self.seq_type]).issuperset(self.seq)

    def get_bioseq_type(self):
        """Returns sequence type"""
        return self.seq_type

    def get_seq_info(self):
        """Returns 4 strings & full sequence info"""
        return f"[Label]: {self.label}\n[Sequence]: {self.seq}\n[Biotype]: {self.seq_type}\n[Length]: {len(self.seq)}"

    def generate_rnd_seq(self, length=10, seq_type="DNA"):
        """Generate a random DNA sequence, provided the length"""
        seq = ''.join([random.choice(NUCLEOTIDE_BASE[seq_type])
                       for x in range(length)])
        self.__init__(seq, seq_type, "Randomly generated sequence")

    def nuc_freq(self):
        """To count the frequency of each nucleotide in the seq & return a dict"""
        return dict(Counter(self.seq))

    def transcription(self):
        """DNA->RNA Transcription. Replacing Thymine with Uracil"""
        if self.seq_type == "DNA":
            return self.seq.replace("T", "U")
        return "Not a DNA sequence"

    def reverse_complement(self):
        """Swapping A with T and G with C. Reversing newly generated string."""
        if self.seq_type == "DNA":
            mapping = str.maketrans("ATCG", "TAGC")
        else:
            mapping = str.maketrans("AUCG", "UAGC")
        return self.seq.translate(mapping)[::-1]

    def gc_content(self):
        """GC Content in a DNA/RNA sequence"""
        return round((self.seq.count("C") + self.seq.count("G")) / len(self.seq) * 100)

    def gc_content_subseq(self, k=20):
        """GC content in a DNA/RNA sub-sequence of length k. k = 15 by default."""
        res = []
        for i in range(0, len(self.seq) - k + 1, k):
            subseq = self.seq[i:i + k]
            res.append(round((subseq.count("C") + subseq.count("G")) / len(subseq) * 100))
        return res

    def translate_seq(self, init_pos=0):
        """Translates DNA sequence into and Amino Acid sequence"""
        if self.seq_type == "DNA":
            return [DNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq) - 2, 3)]
        elif self.seq_type == "RNA":
            return [RNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq) - 2, 3)]

    def codon_freq(self, aminoacid):
        """Provides the frequency of each codon encoding a given amino acid in a DNA sequence"""
        tempList = []
        if self.seq_type == "DNA":
            for i in range(0, len(self.seq) - 2, 3):
                if DNA_Codons[self.seq[i:i + 3]] == aminoacid:
                    tempList.append(self.seq[i:i + 3])
        elif self.seq_type == "RNA":
            for i in range(0, len(self.seq) - 2, 3):
                if RNA_Codons[self.seq[i:i + 3]] == aminoacid:
                    tempList.append(self.seq[i:i + 3])

        freqDict = dict(Counter(tempList))
        totalWeight = sum(freqDict.values())
        for seq in freqDict:
            freqDict[seq] = round(freqDict[seq] / totalWeight, 2)
        return freqDict

    def gen_reading_frames(self):
        frame = []
        frame.append(self.translate_seq(0))
        frame.append(self.translate_seq(1))
        frame.append(self.translate_seq(2))
        tmp_seq = bio_seq(self.reverse_complement(), self.seq_type)
        frame.append(tmp_seq.translate_seq(0))
        frame.append(tmp_seq.translate_seq(1))
        frame.append(tmp_seq.translate_seq(2))
        del tmp_seq
        return frame

    def proteins_from_rf(self, aa_seq):
        """Compute all possible proteins in an amino acid sequence and return a list of possible proteins"""
        current_prot = []
        protein = []
        for aa in aa_seq:
            if aa == "_":
                if current_prot:
                    for p in current_prot:
                        protein.append(p)
                    current_prot = []
            else:
                if aa == "M":
                    current_prot.append("")
                for i in range(len(current_prot)):
                    current_prot[i] += aa
        return protein

    def all_prots_from_orfs(self, startReadPos=0, endReadPos=0, ordered=False):
        """Compute all possible proteins for all open reading frames"""
        """Protein search DB: https://www.ncbi.nlm.nih.gov/nuccore/NM_001185097.2"""
        """API can be used to pull protein info"""
        if endReadPos > startReadPos:
            tmp_seq = bio_seq(
                self.seq[startReadPos: endReadPos], self.seq_type)
            rfs = tmp_seq.gen_reading_frames()
        else:
            rfs = self.gen_reading_frames()

        res = []
        for rf in rfs:
            prots = self.proteins_from_rf(rf)
            for p in prots:
                res.append(p)

        if ordered:
            return sorted(res, key=len, reverse=True)
        return res

