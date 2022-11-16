import requests
import json


class EbiAPI():
    '''
    commute with ebi's api
    https://www.ebi.ac.uk/proteins/api/doc/index.html
    '''

    def __init__(self, uniprot_id):

        #get protein info from ebi
        requestURL = "https://www.ebi.ac.uk/proteins/api/proteins/{}".format(
            uniprot_id)
        r = requests.get(requestURL, headers={"Accept": "application/json"})
        r.raise_for_status()

        #to json
        r = json.loads(r.text)

        #get some info for further use
        self.json_raw = r
        self.uniprot_id = r['accession']
        self.protein_sequence = r['sequence']['sequence']
        self.taxonomy = r['organism']['taxonomy']
        self.gene_name = r['gene'][0]['name']['value']
