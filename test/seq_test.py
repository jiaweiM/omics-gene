from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from pytest import approx


def test_seq():
    my_seq = Seq('GATCG')
    for index, letter in enumerate(my_seq):
        print("%i %s" % (index, letter))
    print(len(my_seq))


def test_count():
    a_count = Seq("AAAA").count("AA")
    assert a_count == 2


def test_g_count():
    my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
    assert len(my_seq) == 32

    g_count = my_seq.count("G")
    assert g_count == 9

    gc_percentage = 100 * float(my_seq.count('G') + my_seq.count('C')) / len(my_seq)
    assert gc_percentage == approx(46.875)

    gc = GC(my_seq)
    assert gc == approx(46.875)


def test_ctr():
    my_seq = Seq("AGTACACTGGT")
    comp_seq = my_seq.complement()
    assert str(comp_seq) == 'TCATGTGACCA'

    rev_comp_seq = my_seq.reverse_complement()
    assert str(rev_comp_seq) == 'ACCAGTGTACT'


def test_slice():
    my_seq = Seq('GATCGATGGGCCTATATAGGATCGAAAATCGC')
    sub_seq = my_seq[4:12]
    assert str(sub_seq) == 'GATGGGCC'


def test_concatenate():
    protein_seq = Seq('EVRNAK')
    dna_seq = Seq('ACGT')
    plus_seq = protein_seq + dna_seq
    assert str(plus_seq) == 'EVRNAKACGT'

    list_of_seqs = [Seq("ACGT"), Seq("AACC"), Seq("GGTT")]
    concatenated_seq = Seq('')
    for s in list_of_seqs:
        concatenated_seq += s
    assert str(concatenated_seq) == 'ACGTAACCGGTT'


def test_join():
    contigs = [Seq("ATG"), Seq("ATCCCG"), Seq("TTGCA")]
    spacer = Seq('N' * 3)
    all_seq = spacer.join(contigs)
    assert str(all_seq) == 'ATGNNNATCCCGNNNTTGCA'


def test_complement():
    my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
    comp_seq = my_seq.complement()
    assert str(comp_seq) == 'CTAGCTACCCGGATATATCCTAGCTTTTAGCG'


def test_transcription():
    coding_dna = Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')
    template_dna = coding_dna.reverse_complement()
    assert str(template_dna) == 'CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT'

    message_rna = coding_dna.transcribe()
    assert str(message_rna) == 'AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG'

    message_rna2 = coding_dna.replace('T', 'U')
    assert message_rna == message_rna2


def test_translate():
    message_rns = Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')
    protein_seq = message_rns.translate()
    print(protein_seq)


def test_translation_table():
    standard_table = CodonTable.unambiguous_dna_by_name['Standard']
    print(standard_table)
    mito_table = CodonTable.unambiguous_dna_by_name['Vertebrate Mitochondrial']
