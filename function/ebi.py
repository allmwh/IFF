import requests
import json

class ProteinNotFound(Exception):
    pass

class EbiAPI():
    '''
    communicate with EBI's API to get sequence infomations
    https://www.ebi.ac.uk/proteins/api/doc/index.html
    '''

    def __init__(self, uniprot_id):

        # get protein info from EBI
        requestURL = "https://www.ebi.ac.uk/proteins/api/proteins/{}".format(
            uniprot_id)
        
        try:
            r = requests.get(requestURL, headers={"Accept": "application/json"})
            r.raise_for_status()
        except:
            raise ProteinNotFound("Please check Uniprot Entry ID")

        # to json
        r = json.loads(r.text)

        # get more information for further use
        self.json_raw = r
        self.uniprot_id = r['accession']
        self.protein_sequence = r['sequence']['sequence']
        self.taxonomy = r['organism']['taxonomy']
        self.gene_name = r['gene'][0]['name']['value']
