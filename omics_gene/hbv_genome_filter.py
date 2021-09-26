from Bio import SeqIO
from Bio.Seq import Seq

gene_list = list(SeqIO.parse(r"D:\data\hbv\hbv_all.gb", "genbank"))

print("Total genome count: ", len(gene_list))

# 只需要环状 DNA，移除非环状的序列
circular_gene_list = []
for record in gene_list:
    annotations = record.annotations
    topology = annotations['topology']
    if topology == 'circular':
        circular_gene_list.append(record)

print("Circular genome count: ", len(circular_gene_list))


# SeqIO.write(circular_gene_list, r'D:\data\hbv\hbv_circular.gb', 'genbank')

def is_b_type(a_record):
    for feature in a_record.features:
        if feature.type == 'source':
            qualifiers = feature.qualifiers
            if 'note' in qualifiers:
                notes = qualifiers['note']
                for note in notes:
                    if 'genotype: B' in note:
                        return True
    return False


def is_c_type(a_record):
    for feature in a_record.features:
        if feature.type == 'source':
            qualifiers = feature.qualifiers
            if 'note' in qualifiers:
                notes = qualifiers['note']
                for note in notes:
                    if 'genotype: C' in note:
                        return True
    return False


def is_type(a_record, name):
    for feature in a_record.features:
        if feature.type == 'source':
            qualifiers = feature.qualifiers
            if 'note' in qualifiers:
                notes = qualifiers['note']
                for note in notes:
                    if name in note:
                        return True
    return False


type_a_list = []
type_b_list = []
type_c_list = []
type_d_list = []
type_e_list = []
type_f_list = []
type_g_list = []
type_h_list = []
type_i_list = []

for record in circular_gene_list:
    if is_type(record, 'genotype: A'):
        type_a_list.append(record)
    elif is_type(record, 'genotype: B'):
        type_b_list.append(record)
    elif is_type(record, "genotype: C"):
        type_c_list.append(record)
    elif is_type(record, 'genotype: D'):
        type_d_list.append(record)
    elif is_type(record, 'genotype: E'):
        type_e_list.append(record)
    elif is_type(record, 'genotype: F'):
        type_f_list.append(record)
    elif is_type(record, 'genotype: G'):
        type_g_list.append(record)
    elif is_type(record, 'genotype: H'):
        type_h_list.append(record)
    elif is_type(record, 'genotype: I'):
        type_i_list.append(record)

print("HBV genotype a count: ", len(type_a_list))
print("HBV genotype b count: ", len(type_b_list))
print("HBV genotype c count: ", len(type_c_list))
print("HBV genotype d count: ", len(type_d_list))
print("HBV genotype e count: ", len(type_e_list))
print("HBV genotype f count: ", len(type_f_list))
print("HBV genotype g count: ", len(type_g_list))
print("HBV genotype h count: ", len(type_h_list))
print("HBV genotype i count: ", len(type_i_list))

SeqIO.write(type_a_list, r'D:\data\hbv\hbv_genptype_a.fa', 'fasta')
SeqIO.write(type_b_list, r'D:\data\hbv\hbv_genptype_b.fa', 'fasta')
SeqIO.write(type_c_list, r'D:\data\hbv\hbv_genptype_c.fa', 'fasta')
SeqIO.write(type_d_list, r'D:\data\hbv\hbv_genptype_d.fa', 'fasta')
SeqIO.write(type_e_list, r'D:\data\hbv\hbv_genptype_e.fa', 'fasta')
SeqIO.write(type_f_list, r'D:\data\hbv\hbv_genptype_f.fa', 'fasta')
SeqIO.write(type_g_list, r'D:\data\hbv\hbv_genptype_g.fa', 'fasta')
SeqIO.write(type_h_list, r'D:\data\hbv\hbv_genptype_h.fa', 'fasta')
SeqIO.write(type_i_list, r'D:\data\hbv\hbv_genptype_i.fa', 'fasta')
