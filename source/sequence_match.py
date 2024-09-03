import numpy as np
import pandas as pd


from Bio import pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def find_best_alignment(query_sequence: str, reference_sequences: list) -> tuple:
    """
    Aligns a query sequence to a list of reference sequences and returns the name of the best-aligned sequence.
    
    Parameters:
    query_sequence (str): The sequence of interest that needs to be aligned.
    reference_sequences (list): A list of SeqRecord objects from Biopython, containing sequences and their names.
    
    Returns:
    tuple: The name of the most aligned reference sequence and best alignment object
    """
    best_score = -1
    best_match_name = None
    best_alignment = None


    for record in reference_sequences:
        alignments = pairwise2.align.globalxx(query_sequence, record.seq)
        # Take the alignment with the highest score
        score = alignments[0][2]  # [2] contains the alignment score
        
        if score > best_score:
            
            best_score = score
            best_match_name = record.id
            best_alignment = alignments[0]

    
    return best_match_name, best_alignment



