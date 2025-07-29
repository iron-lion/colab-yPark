import requests
from collections import defaultdict, Counter
from Bio import SeqIO # You'll need Biopython installed (pip install biopython)
from Bio import motifs
from Bio.Seq import Seq
import tqdm


class JASPARrest:
    def __init__(self, use_progress_bar=False):
        self.endpoint = 'https://jaspar.elixir.no'
        self.headers = {"content-type": "application/json"}
        self.use_progress_bar = use_progress_bar

    def get_pfm(self, taxid):
        _jaspar_motifs = dict()

        url = f"{self.endpoint}/api/v1/taxon/{taxid}?page_size=1205"
        response = requests.get(url, headers=self.headers)
        if response.ok:
            data = response.json()['results']
            matrix_ids = [x['matrix_id'] for x in data]

            pfm_url = f"{self.endpoint}/api/v1/matrix"
            for mid in tqdm.tqdm(matrix_ids):
                r = requests.get(f"{pfm_url}/{mid}", headers=self.headers)
                pfm_matrix = r.json()
                if r.ok:
                    current_motif_name = f"{pfm_matrix['matrix_id']},{pfm_matrix['name']}"
                    motif_obj = motifs.Motif(alphabet=["A","C","G","T"], counts=pfm_matrix['pfm'])

                    motif_obj.name = current_motif_name
                    _jaspar_motifs[current_motif_name] = motif_obj

            return _jaspar_motifs
        else:
            raise('Connection to JASPAR, failed!')

    def browse(self, taxid):
        url = f"{self.endpoint}/api/v1/tffm/{taxid}"
        response = requests.get(url, headers=self.headers)
        print(response.json())


# --- Main Execution ---
if __name__ == "__main__":
    jaspar_taxgrp = {"C.elegans": "Nematodes",
                    "Killifish": "Vertebrates",
                    "Mouse": "Vertebrates",
    }
    target_taxid = jaspar_taxgrp['C.elegans']

    rest_jaspar = JASPARrest()
    jaspar_motifs = rest_jaspar.get_pfm(target_taxid)
    print(jaspar_motifs)
