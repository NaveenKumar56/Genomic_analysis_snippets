import os
from Bio import SeqIO
import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import shutil

def merge_fasta(input_fasta, output_fasta, gap_size=100):
    seq_strings = []
    records = list(SeqIO.parse(input_fasta, "fasta"))
    gap_string = "N" * gap_size
    for r in records:
        seq_str = str(r.seq)
        seq_strings.append(seq_str)
    merged_seq = gap_string.join(seq_strings)
    merged_record = SeqRecord(merged_seq, id="merged_contig",
                              description=f"All contigs merged with {gap_size} N gap")
    SeqIO.write(merged_record, output_fasta, "fasta")



def rename_fasta_header(output_fasta, output_renamed_fasta):
    file_name = os.path.splitext(os.path.basename(output_fasta))[0]
    records = list(SeqIO.parse(output_fasta, "fasta"))

    for r in records:
        r.id = file_name
        r.name = file_name
        r.description = ""

    # Ensure output folder exists
    SeqIO.write(records, output_renamed_fasta, "fasta")


def main(input_folder, working_directory):
    path = input_folder
    complete_genome_path = os.path.join(working_directory, "complete_genomes")
    output_renamed_fasta =  os.path.join(working_directory, "complete_renamed_genomes")
    os.makedirs(complete_genome_path, exist_ok=True)
    os.makedirs(output_renamed_fasta, exist_ok=True)

    single_copy_genome = 0
    multi_copy_genome = 0
    for file in os.listdir(path):
    

        file_path = os.path.join(path, file)
    
        if not os.path.isfile(file_path):
            continue  # skip directories
    
        records = list(SeqIO.parse(file_path, "fasta"))
        num_contigs = len(records)
    
    # Use the original filename as is, without changing extension
        output_file = os.path.join(complete_genome_path, file)
    
        if num_contigs == 1:
            shutil.copy2(file_path, output_file)
            single_copy_genome += 1
            print(f"Copied single-contig file: {file}")
        else:
            multi_copy_genome += 1
            merge_fasta(file_path, output_file, gap_size=100)
            print(f"Merged multi-contig file: {file}")

    print(f"File: {file} | Single-copy genomes: {single_copy_genome} | Multi-copy genomes: {multi_copy_genome}")



    for file in os.listdir(complete_genome_path):
        file_path = os.path.join(complete_genome_path, file)
        output_file = os.path.join(output_renamed_fasta, file)
        rename_fasta_header(file_path, output_file)
        print(f"Renamed headers in file: {file}")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge FASTA contigs and rename headers to filename")
    parser.add_argument("input_folder", help="Path to input folder containing FASTA files")
    parser.add_argument("working_directory", help="Path to folder where merged and renamed FASTA files will be stored")

    args = parser.parse_args()
    main(args.input_folder, args.working_directory)
