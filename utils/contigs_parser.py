from pathlib import Path
from motif_utils import seq2kmer
import pandas as pd
import os


def create_full_path(processed_dir, extension=".tsv"):
    return (processed_dir/Path("dev")).with_suffix(extension)

def get_phages_labels(size):
    return size*[0]

def get_bacteria_labels(size):
    return size*[1]

if __name__ == "__main__":
    """
    Splits contig.fasta file to nodes, each node saves at different .txt file
    """
    K_MER_SIZE = 6
    data_dir = Path("/home/danya-sakharov/dna_bert/data/bacterias/one")
    for file in os.listdir(data_dir):
        if ".fasta" in file:
            filepath = data_dir / file
            with open(filepath) as f:
                lines = [line.rstrip() for line in f]

            processed_dir = data_dir / filepath.stem / f"K{K_MER_SIZE}"

            if not processed_dir.exists():
                processed_dir.mkdir(parents=True, exist_ok=True)

            full_contigs = []
            contig = []
            label = []
            node_infos = [lines[0]]

            for line in lines[1:]:
                line = line.upper()
                if ">" in line:
                    contig_joined = "".join(contig)
                    k_mer_contig = seq2kmer(contig_joined, K_MER_SIZE)
                    full_contigs.append(k_mer_contig)
                    contig = []
                    node_infos.append(line)
                else:
                    contig.append(line)
    
            contig_joined = "".join(contig)
            k_mer_contig = seq2kmer(contig_joined, K_MER_SIZE)
            full_contigs.append(k_mer_contig)

            label = get_phages_labels(len(full_contigs))
            print(F"contigs found: {len(full_contigs)}")
            dataset_45 = pd.DataFrame({'sequence':full_contigs, 'label':label, 'node':node_infos})
            dataset_45.to_csv(create_full_path(processed_dir), index=False, sep='\t')
            
            
