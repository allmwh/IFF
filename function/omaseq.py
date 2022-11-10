import re
import json
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path
from tqdm.notebook import trange
from Bio.SeqRecord import SeqRecord

from function.utilities import fasta_to_seqlist

class FetchOmaSeqBatch():
    '''
    get homologs by OMA Group from OMA to fasta, e.g.
    1. get OMA Group fasta from https://omabrowser.org/oma/omagroup/859990/fasta/ 
    2. change sequence name's format, infos are from https://omabrowser.org/api/group/859990/ 
    '''

    def __init__(self):
        pass

    def get_oma_seq(self, oma_group_id, path):
        '''
        pipeline: get fasta from OMA, changing sequence name's format

        oma_group_id: str, OMA Group id, e.g. 859990
        path: str, path for saving fasta

        retrun: None
        '''

        oma_path = Path(path)
        oma_fasta_path = oma_path / "{}.fasta".format(oma_group_id)

        # get raw fasta
        self.__get_oma_fasta(oma_group_id, oma_fasta_path)

        # get fasta info
        fasta_info_dict = self.__get_fasta_info(oma_group_id)

        # get mod info's fasta
        self.__mod_fasta_info(oma_fasta_path, oma_fasta_path, fasta_info_dict)

    def __get_oma_fasta(self, oma_group_id, fasta_path):
        '''
        get OMA Group to fasta file

        oma_group_id: str, OMA Group id, e.g. 859990
        fasta_path: str, path for saving fasta

        return: None
        '''

        try:
            url = "https://omabrowser.org/oma/omagroup/{}/fasta/".format(oma_group_id)
            resp = requests.get(url)
            resp.raise_for_status()
            with open(fasta_path, "w") as file:
                file.write(resp.text)
        except:
            raise Exception("{} get fasta failed from OMA".format(oma_group_id))

    def __get_fasta_info(self, oma_group_id):
        '''
        get sequence infos for each sequence in OMA Group

        oma_group_id: str, OMA Group id, e.g. 859990

        return: dict, fasta_info_dict, contain all sequence's info
        '''

        try:
            url = "https://omabrowser.org/api/group/{}/".format(oma_group_id)
            resp = requests.get(url)
            resp.raise_for_status()
            oma_raw = json.loads(resp.text)
            fasta_info_dict = {}
            for i in oma_raw['members']:

                species = i["species"]["species"]
                
                #sometimes species name are too long, remove some infos like strain
                species = re.sub("\(.*\)", "", species) 
                oma_id = i["omaid"]
                
                canonical_id = i["canonicalid"]
                taxon_id = i["species"]["taxon_id"]

                fasta_info_dict[oma_id] = {
                    "oma_id": oma_id,
                    "species": species,
                    "canonical_id": canonical_id,
                    "taxon_id": taxon_id,
                }

            return fasta_info_dict
        except:
            raise Exception("{} OMA fetch fasta seqeuence info failed".format(oma_group_id))

    def __mod_fasta_info(self, oma_fasta_path, mod_fasta_path, fasta_info_dict):
        '''
        change sequence name's format

        oma_fasta_path: str, path for saved fasta from __get_oma_fasta()
        mod_fasta_path: str, path for saving fasta 
        fasta_info_dict: dict, fasta sequence's info from __get_fasta_info()

        return: None
        '''

        fasta_list = list(SeqIO.parse(str(oma_fasta_path), 'fasta'))
        mod_fasta_list = []
        for seq_record in fasta_list:
            id = seq_record.id
            record = SeqRecord(seq=seq_record.seq,
                               id=id,
                               description="| {} | {} | {}".format(fasta_info_dict[id]["species"],
                                                                   fasta_info_dict[id]["taxon_id"],
                                                                   fasta_info_dict[id]["canonical_id"])
                            )
            mod_fasta_list.append(record)
        SeqIO.write(mod_fasta_list, mod_fasta_path, "fasta")

class TaxSeqFilter():
    '''
    filter homologs to newer fasta by taxonomy id
    '''

    def __init__(self, taxonomy):
        '''
        taxonomy: int, taxonomy id from NCBI for filter
                       NCBI: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=info&id=9606
        '''

        resp = requests.get("https://omabrowser.org/api/taxonomy/{}".format(taxonomy))
        self.taxonomy = taxonomy
        self.taxonomy_list = resp.text

    def taxfilter(self, oma_fasta_path, grouped_fasta_path):
        '''
        oma_fasta_path: str, fasta file path for all OMA paralogs
        grouped_fasta_path: str, fasta file path for grouped paralogs
        
        return: None
        '''

        # read
        oma_fasta_list = fasta_to_seqlist(oma_fasta_path)

        # filter
        filtered_list = []
        for i in oma_fasta_list:
            tax_id = i.description.split("|")[2].replace(" ", "")
            if tax_id in self.taxonomy_list:
                filtered_list.append(i)

        with open(grouped_fasta_path, "w") as output_handle:
            SeqIO.write(filtered_list, output_handle, "fasta")