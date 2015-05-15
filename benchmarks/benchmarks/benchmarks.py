# Write the benchmarking functions here.
# See "Writing benchmarks" in the asv docs for more information.
import skbio
from qiime_default_reference import (get_template_alignment,
                                    get_reference_sequences)


class TimeSuite:
    """
    Benchmarks for the sequence class
    """
    def setup(self):
        # load data from a file
        self.gapped_sequences = [(s.id, str(s)) for s in skbio.SequenceCollection.read(get_template_alignment())][:500]
        self.sequences = [(s.id, str(s)) for s in skbio.SequenceCollection.read(get_reference_sequences())][:500]
        self.motif_1 = "GGTGCAAGCCGGTGGAAACA"
        self.DNA_sequences = [skbio.sequence.DNA(seq, id=str(id)) for id_, seq in self.sequences]
        self.gapped_DNA_sequences = [skbio.sequence.DNA(seq, id=str(id)) for id_, seq in self.gapped_sequences]
        self.sgc = skbio.sequence.genetic_code(1)
        self.motif = "(GGTGCAAGCCGGTGGAAACA)"
        self.empty_list = []

    def time_reverse_complement(self):
        for seq in self.DNA_sequences:
            seq.reverse_complement()

    def time_object_creation(self):
        for id_, seq in self.sequences:
            skbio.sequence.DNA(seq, id=id_, validate=False)

    def time_object_creation_validate(self):
        for id_, seq in self.sequences:
            skbio.sequence.DNA(seq, id=id_)

    def time_degap_all(self):
        for seq in self.gapped_DNA_sequences:
            seq.degap()

    def time_translate(self):
        for seq in self.DNA_sequences:
            self.sgc.translate(seq, 1)

    def time_search_for_motif(self):
        for seq in self.DNA_sequences:
            list(seq.slices_from_regex(self.motif))

    def time_kmer_count_5(self):
        for seq in self.DNA_sequences:
            seq.kmer_frequencies(5)

    def time_kmer_count_25(self):
        for seq in self.DNA_sequences:
            seq.kmer_frequencies(5)

    def time_validate_chars(self):
        for seq in self.DNA_sequences:
            skbio.sequence.DNA(seq)

    def time_filter_invalid_seqs(self):
        for seq in self.DNA_sequences:
            try:
                self.empty_list.append(skbio.sequence.DNA(seq))
            except ValueError:
                pass

    def time_expand_degenerates(self):
        for seq in self.DNA_sequences:
            list(seq.expand_degenerates())

    def time_gc_content(self):
        for seq in self.DNA_sequences:
            float(seq.count("G") + seq.count("C"))/len(seq)

    def time_find_motif_in_gapped(self):
        for seq in self.gapped_DNA_sequences:
            list(seq.slices_from_regex(self.motif, ignore=seq.gaps()))


# class MemSuite:
#     def mem_list(self):
#         return [0] * 256
