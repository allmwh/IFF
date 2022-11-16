from pathlib import Path
from tqdm.notebook import trange
from Bio.Align.Applications import ClustalOmegaCommandline


class Alignment():
    """
    alignment by Biopython ClustalOmegaCommandline
    plese install ClustalOmega: http://www.clustal.org/omega/
    Ubuntu/Debian based linux distro: sudo apt install clustalo
    """

    def __init__(self):
        pass

    def alignment_single(self, input_fasta, output_fasta, delete=False):
        '''
        fasta alignment

        input_fasta: str, fasta file for alignment
        output_fasta: str, fasta file for alied fasta

        return: None
        '''
        input_fasta = Path(input_fasta)
        output_fasta = Path(output_fasta)

        try:
            clustalomega_cline = ClustalOmegaCommandline(infile=input_fasta,
                                                         outfile=output_fasta,
                                                         verbose=True,
                                                         auto=True,
                                                         force=True)
            clustalomega_cline()
            # remove input file if successfully alignment
            if delete:
                input_fasta.unlink()

        except:
            raise Exception("{} alignment failed".format(str(input_fasta)))

    def alignment_path(self, input_path, output_path, delete=False):
        '''
        alignment all file under input_path, and save alignment fasta to output_path

        input_path: str, path saving many fasta to be alied 
        output_path: str, path for alied fsta output

        return: None
        '''

        failed_list = []

        input_path = Path(input_path)
        input_path.mkdir(parents=True, exist_ok=True)
        output_path = Path(output_path)
        output_path.mkdir(parents=True, exist_ok=True)

        fasta_pathlist = list(Path(input_path).rglob("*.fasta"))

        t = trange(len(fasta_pathlist), desc="", leave=True, position=0)
        for i in t:
            fasta_path = fasta_pathlist[i].parts[-1]
            t.set_description(fasta_path)

            input_fasta = input_path / fasta_path
            output_fasta = output_path / fasta_path

            try:
                self.alignment_single(input_fasta, output_fasta, delete=delete)
            except Exception as e:
                print(e)
                failed_list.append(fasta_path)

        return failed_list
