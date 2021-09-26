from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def test_ctr():
    simple_seq = Seq("GATC")
    simple_record = SeqRecord(simple_seq)
    assert simple_record.id == '<unknown id>'
    assert simple_record.description == '<unknown description>'
    assert simple_record.seq == Seq('GATC')
