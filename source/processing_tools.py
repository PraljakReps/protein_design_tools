import pandas as pd
from Bio import SeqIO

def fasta_to_csv(fasta_file: str, csv_file: str, header_col: str = "Header", seq_col: str = "Sequence"):
    """
    Converts a FASTA file to a CSV file.

    Parameters:
    fasta_file (str): Path to the input FASTA file.
    csv_file (str): Path to the output CSV file.
    header_col (str): Name of the column for sequence headers in the CSV file.
    seq_col (str): Name of the column for sequences in the CSV file.
    """
    records = list(SeqIO.parse(fasta_file, "fasta"))
    data = {
        header_col: [record.id for record in records],
        seq_col: [str(record.seq) for record in records]
    }
    df = pd.DataFrame(data)
    df.to_csv(csv_file, index=False)
    print(f"FASTA file has been successfully converted to CSV and saved as {csv_file}")

def csv_to_fasta(csv_file: str, fasta_file: str, header_col: str = "Header", seq_col: str = "Sequence"):
    """
    Converts a CSV file to a FASTA file.

    Parameters:
    csv_file (str): Path to the input CSV file.
    fasta_file (str): Path to the output FASTA file.
    header_col (str): Name of the column containing sequence headers in the CSV file.
    seq_col (str): Name of the column containing sequences in the CSV file.
    """
    df = pd.read_csv(csv_file)
    with open(fasta_file, "w") as output_handle:
        for _, row in df.iterrows():
            output_handle.write(f">{row[header_col]}\n{row[seq_col]}\n")
    print(f"CSV file has been successfully converted to FASTA and saved as {fasta_file}")
