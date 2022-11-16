import re
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq


def seq_aa_check(sequence):
    '''
    check 20 amino acid abbreviations are usual used and delete unvalied letter 
    ref: https://www.ncbi.nlm.nih.gov/Class/MLACourse/Modules/MolBioReview/iupac_aa_abbreviations.html
    
    sequence: str, seqeunce

    return: str, sequence with checked aa
    '''

    sequence = re.sub('[^ARNDCEQGHILKMFPSTWYV]', '', sequence)

    return sequence


def fasta_seq_aa_check(fasta_path, checked_fasta_path):
    '''
    check 20 amino acid within fasta file

    fasta_path: str, fasta file to check
    checked_fasta_path: str, checked fasta file

    return:None
    '''

    fasta_path = Path(fasta_path)
    fasta_path = fasta_path.absolute()
    checked_fasta_path = Path(checked_fasta_path)
    checked_fasta_path = checked_fasta_path.absolute()

    # read former fasta
    fasta_list = list(SeqIO.parse(fasta_path, "fasta"))

    # seq aa check
    for index, element in enumerate(fasta_list):
        former_sequence = str(element.seq)
        aa_checked_sequence = Seq(seq_aa_check(former_sequence))
        fasta_list[index].seq = aa_checked_sequence

    # write to new path
    _ = SeqIO.write(fasta_list, checked_fasta_path, "fasta")


def get_fasta_filename(path):
    '''
    simple function to get fasta file name
    '''
    path = Path(path)
    filename = path.parts[-1].split(".")[0]
    return filename


def fasta_to_seqlist(path):
    '''
    simple function for loading fasta file to list of BioPython sequence object

    path: path for fasta file

    return: list with BioPython sequence object
    '''
    return list(SeqIO.parse(str(path), 'fasta'))


def get_fasta_seq_info(path):
    '''
    get fasta sequence info, including oma_group_id(OMA groups id from fasta file name), 
    all tax sequences info in fasta
    '''
    path = Path(path)
    path = path.absolute()
    fasta_list = list(SeqIO.parse(path, "fasta"))
    seq_info_list = []
    for i in fasta_list:
        sequence_info = i.description
        sequence_info = sequence_info.split('|')

        oma_protein_id = sequence_info[0].strip()
        species = sequence_info[1].strip()
        taxon_id = sequence_info[2].strip()
        oma_cross_reference = sequence_info[3].strip()

        seq_info_list.append({
            "oma_protein_id": oma_protein_id,
            "species": species,
            "taxon_id": taxon_id,
            "oma_cross_reference": oma_cross_reference
        })

    oma_group_id = path.parts[-1].split(".")[0]

    return {"oma_group_id": oma_group_id, "homologous_info": seq_info_list}


def find_human_sequence(path):
    """
    find human sequence from fasta file
    """
    path = Path(path)
    all_sequence_list = fasta_to_seqlist(path)
    oma_group_id = get_fasta_seq_info(path)['oma_group_id']

    for i in all_sequence_list:
        tax_id = i.description.split("|")[2]
        if int(tax_id) == 9606:
            sequence = str(i.seq)
            return {
                "oma_group_id": oma_group_id,
                "sequence_name": i.description,
                "remove_gap_sequence": sequence.replace("-", ""),
                "sequence": sequence
            }
