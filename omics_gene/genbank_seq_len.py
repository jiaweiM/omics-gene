from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def reorder_seq(seq, order_len):
    return seq[order_len:] + seq[:order_len]


def test_reorder():
    record = SeqRecord(Seq("AACTCCACAAGTGG"), id='AY220704.1')
    new_record = reorder_seq(record, 2)
    assert new_record.seq == Seq('CTCCACAAGTGGAA')


def test_length(record_list, length=3215):
    retain_record_list = []
    for record in record_list:
        if len(record) == length:
            retain_record_list.append(record)
    return retain_record_list


def stat_length(record_list):
    record_length_dict = defaultdict(int)
    for record in record_list:
        record_length_dict[len(record)] += 1
    return record_length_dict


def stat_n(record_list):
    """
    统计序列中 N 的个数
    :param record_list:
    :return:
    """
    n_dict = defaultdict(int)
    for record in record_list:
        n_count = record.seq.count('N')
        n_dict[n_count] += 1
    return n_dict


def stat_start(record_list):
    n_dict = defaultdict(int)
    for record in record_list:
        start_seq = str(record.seq[:9])
        n_dict[start_seq] += 1
    return n_dict


def stat_end(record_list):
    n_dict = defaultdict(int)
    for record in record_list:
        end_seq = str(record.seq[-9:])
        n_dict[end_seq] += 1
    return n_dict


def process_genotype_b():
    file = r'D:\data\hbv\hbv_genptype_b.gb'
    length = 3215

    record_list = list(SeqIO.parse(file, 'genbank'))
    print(f'Record in file: {len(record_list)}')

    len_dict = stat_length(record_list)
    print("Record length distribution")
    for record_len in len_dict.keys():
        print(f"{record_len} = {len_dict[record_len]}")

    record_list = test_length(record_list)
    print(f'Record with length {length}: {len(record_list)}')

    n_dict = stat_n(record_list)
    print("Number of N distribution")
    for count in n_dict.keys():
        print(f'{count} = {n_dict[count]}')

    start_dict = stat_start(record_list)
    print("Start seq distribution")
    for sub_seq in start_dict.keys():
        print(f"{sub_seq} = {start_dict[sub_seq]}")

    end_dict = stat_end(record_list)
    print('End seq distribution')
    for sub_seq in end_dict.keys():
        print(f"{sub_seq} = {end_dict[sub_seq]}")

    # reorder_record = []
    # with_n_record = 0
    # aa_record = 0
    # for record in record_list:
    #     if record.seq.count('N') > 0:
    #         with_n_record += 1
    #         continue
    #
    #     # CTCCACCAC
    #     if record.seq.startswith('CTCCACCAC'):
    #         reorder_record.append(record)
    #         continue
    #
    #     if record.seq.startswith('AA'):
    #         aa_record += 1
    #         reorder_record.append(reorder_seq(record, 2))
    #         continue
    #
    #     index = record.find('CTCCACCAC')


process_genotype_b()

# SeqIO.write(retain_record_list, 'D:\data\hbv\hbv_genptype_b_len.gb', 'fasta')
