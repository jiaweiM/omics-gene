from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.Seq import Seq
import pytest
from pytest import approx


@pytest.fixture()
def record():
    record = next(SeqIO.parse('NC_003977.gb', 'genbank'))
    return record


def test_id(record):
    genbank_id = record.id
    assert genbank_id == 'NC_003977.2'  # VERSION 部分


def test_name(record):
    gb_name = record.name  # LOCUS id
    assert gb_name == 'NC_003977'


def test_description(record):
    assert record.description == 'Hepatitis B virus (strain ayw) genome'  # DEFINITION


def test_letter_annotation(record):
    letter_annotations = record.letter_annotations  # empty
    assert not letter_annotations


def test_annotations(record):
    annotations = record.annotations
    assert annotations['molecule_type'] == 'DNA'  # mol_type
    assert annotations['topology'] == 'circular'
    assert annotations['data_file_division'] == 'VRL'
    assert annotations['date'] == '25-MAR-2021'
    assert annotations['accessions'] == ['NC_003977']
    assert annotations['sequence_version'] == 2
    assert annotations['keywords'] == ['RefSeq', 'antigen', 'core antigen', 'genome', 'signal peptide']
    assert annotations['source'] == 'Hepatitis B virus'
    assert annotations['organism'] == 'Hepatitis B virus'
    assert annotations['taxonomy'] == ['Viruses', 'Riboviria', 'Pararnavirae', 'Artverviricota', 'Revtraviricetes',
                                       'Blubervirales', 'Hepadnaviridae', 'Orthohepadnavirus']
    assert annotations[
               'comment'] == 'VALIDATED REFSEQ: This record has undergone validation or\npreliminary review. The reference sequence is identical to V01460.\nOn Oct 23, 2015 this sequence version replaced NC_003977.1.\nCOMPLETENESS: full length.'


def test_features(record):
    features = record.features
    assert len(features) == 28
    feature1 = features[0]
    assert feature1.type == 'source'
    print(feature1)
    location = feature1.location



def test_gc():
    my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
    ratio = GC(my_seq)
    assert ratio == approx(46.875)


def test_seq(record):
    seq = record.seq
    assert len(seq) == 3182
    assert seq[0] == 'A'

    a_count = seq.count('A')
    assert a_count == 731

    print(a_count)
