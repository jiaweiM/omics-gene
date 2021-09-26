from Bio.Seq import Seq
import tensorflow as tf


def test_seq():
    dna = Seq('CTCCCTTATCGTCAATCTTCTCG')
    print(dna.reverse_complement_rna())
    print(dna.transcribe())
