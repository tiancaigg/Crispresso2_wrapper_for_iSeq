#!/usr/bin/env python
# coding: utf-8

import argparse
import subprocess
import glob
import os
import pandas as pd
from multiprocessing import Pool
from Bio import SeqIO
import traceback
import sys
import logging

# Create the parser and add arguments
parser = argparse.ArgumentParser(description='Run CRISPResso analysis.')
parser.add_argument('--fastq_dir', required=True, help='Directory containing FASTQ files')
parser.add_argument('--amplicon_fasta_dir', required=True, help='Directory containing amplicon FASTA files')
parser.add_argument('--output_dir', default='./Crispresso_output', help='Output directory for CRISPResso (default: ./Crispresso_output)')
parser.add_argument('--amplicon_book', required=True, help='Amplicon booking list in xlsx format')
parser.add_argument('--sheet', required=True, help='Sheet name in the Amplicon booking list')
parser.add_argument('--verbose', type=str, choices=['T', 'F'], default='T', help='Verbose output (T/F, default: T)')

# Parse the arguments
args = parser.parse_args()

# Use the arguments in the script
fastq_dir = args.fastq_dir
amplicon_fasta_dir = args.amplicon_fasta_dir
output_dir = args.output_dir
amplicon_book_path = args.amplicon_book
sheet_name = args.sheet
verbose = args.verbose == 'T'

if not os.path.exists(output_dir):
    os.mkdir(output_dir) 
else:
    print('Output directory already exists, use --output_dir to change the output directory')
    sys.exit()

# Set up logging
log_file = os.path.join(output_dir, 'crispresso_wrapper.log')
logging.basicConfig(level=logging.INFO if verbose else logging.ERROR,
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    handlers=[
                        logging.FileHandler(log_file),
                        logging.StreamHandler()
                    ])

# Function to log messages both to file and console
def log_message(message, level=logging.INFO):
    if verbose:
        logging.log(level, message)

log_message(f"Starting CRISPResso analysis with arguments: {args}")

amplicon_book = pd.read_excel(amplicon_book_path, sheet_name=sheet_name)

amplicon_book_long = pd.wide_to_long(amplicon_book, stubnames=['user', 'accession/gene', 'PCR_Product_fa', 'Guide'], 
                                     i=['fastq', 'Index1', 'Index2'],
                                     j=' Set', suffix=r' \Set\w')

# Remove all whitespaces 
df_obj = amplicon_book_long.select_dtypes(['object'])
amplicon_book_long[df_obj.columns] = df_obj.apply(lambda x: x.str.strip())
amplicon_book_long = amplicon_book_long.reset_index()

amplicon_book_long = amplicon_book_long[['fastq', 'PCR_Product_fa', 'Guide', 'accession/gene']] 
amplicon_book_long = amplicon_book_long[amplicon_book_long.PCR_Product_fa.notnull()]
amplicon_book_long['PCR_Product_fa'] = amplicon_book_long['PCR_Product_fa'].str.replace('.fa$', '', regex=True)
amplicon_book_long['Guide'] = amplicon_book_long['Guide'].astype(str)
amplicon_book_long['Guide'] = amplicon_book_long['Guide'].str.upper()

amplicon_book_long['Guide'] = amplicon_book_long['Guide'].str.upper().str.strip()

actual_files = set(os.path.basename(f) for f in glob.glob(amplicon_fasta_dir + '*.fa'))
valid_files_mask = amplicon_book_long['PCR_Product_fa'].apply(lambda x: x + '.fa' in actual_files)
invalid_rows = amplicon_book_long[~valid_files_mask]
amplicon_book_long = amplicon_book_long[valid_files_mask]

if not invalid_rows.empty:
    log_message(f'These fasta files are missing: {set(invalid_rows["PCR_Product_fa"])}', logging.WARNING)

def get_amplicon_seq(file_path):
    try:
        seq_record = SeqIO.read(file_path, "fasta")
        return str(seq_record.seq)
    except ValueError as e:
        raise ValueError(f"Error reading FASTA file {file_path}: {str(e)}")

def get_prefix_fastq(s_value):
    number = int(s_value.strip("S"))
    return f"{s_value}_{s_value}" if number >= 10 else f"{s_value}_S{number}"

def run_crispresso(row):
    fq_file, amplicon_fasta, guide = row['fastq'], row['PCR_Product_fa'], row['Guide']
    amplicon_file = os.path.join(amplicon_fasta_dir, amplicon_fasta + '.fa')
    fastq_prefix = get_prefix_fastq(fq_file)
    log_file = os.path.join(output_dir, amplicon_fasta, fastq_prefix + '_' + amplicon_fasta + '_error_log.txt')
    
    log_message(f"Processing sample: {amplicon_fasta}")
    
    try:
        amplicon_seq = get_amplicon_seq(amplicon_file)
        
        crispresso_command = [
            'CRISPResso',
            '--fastq_r1', os.path.join(fastq_dir, fastq_prefix + '_L001_R1_001.fastq.gz'),
            '--fastq_r2', os.path.join(fastq_dir, fastq_prefix + '_L001_R2_001.fastq.gz'),
            '--amplicon_seq', amplicon_seq,
            '--keep_intermediate',
            '--amplicon_name', amplicon_fasta,
            '--bam_output',
            '-o', os.path.join(output_dir, amplicon_fasta),
            '--file_prefix', fastq_prefix + '_' + amplicon_fasta,
        ]
        
        if not pd.isna(guide):
            crispresso_command.extend(['--guide_seq', guide])
        
        log_message(f"Running CRISPResso command for {amplicon_fasta}")
        crispresso_process = subprocess.run(crispresso_command, capture_output=True, text=True)
        if crispresso_process.returncode == 0:
            log_message(f"CRISPResso command was successful for {amplicon_fasta}")
            return f"CRISPResso command was successful for {amplicon_fasta}"
        else:
            error_message = f"CRISPResso command failed for {amplicon_fasta}! Error: {crispresso_process.stderr}"
            log_message(error_message, logging.ERROR)
            with open(log_file, 'a') as f:
                f.write(error_message + "\n")
            return error_message
    except Exception as e:
        error_message = f"An error occurred for {amplicon_fasta}: {str(e)}"
        log_message(error_message, logging.ERROR)
        with open(log_file, 'a') as f:
            f.write(error_message + "\n" + traceback.format_exc() + "\n")
        return error_message

log_message("Starting to process samples")
rows_as_dicts = amplicon_book_long.to_dict('records')
pool = Pool(processes=40) 
results = pool.map(run_crispresso, rows_as_dicts)

# Close the pool and wait for all processes to finish
pool.close()
pool.join()

log_message("Finished processing all samples")

# Print results
for result in results:
    log_message(result)

# The pattern to match files deep in the directory structure.
editing_freq_pattern = '**/*quantification_of_editing_frequency.txt'

# Full pattern path combining base directory and the file pattern.
full_editing_freq_pattern = os.path.join(output_dir, editing_freq_pattern)

# Perform the recursive search to get a list of all matching file paths.
editing_file_paths = glob.glob(full_editing_freq_pattern, recursive=True)

# Create a container for the DataFrames
dfs = []
for file_path in editing_file_paths:
    try:
        # Read the current file into a DataFrame
        current_df = pd.read_csv(file_path, sep='\t')
        
        # Add an identifier for each file
        current_df['Source_File'] = os.path.basename(file_path)
        
        # Append the current DataFrame to the container
        dfs.append(current_df)
    except Exception as e:
        log_message(f"Could not read file {file_path}. Error: {e}", logging.ERROR)
        
if dfs:
    final_df = pd.concat(dfs, ignore_index=True) 
    # Regular expression pattern to match strings like "S01", "S02", ..., "S95", etc.
    pattern = r'(?<!\d)(S\d{2})'

    # Extract the pattern from the 'Source_File' column and create a new column 'Fastq' with the result
    final_df['Fastq'] = final_df['Source_File'].str.extract(pattern, expand=False)
    final_df['Insertions%'] = 100 * final_df['Insertions'] / final_df['Reads_aligned']
    final_df['Deletions%'] = 100 * final_df['Deletions'] / final_df['Reads_aligned']
    final_df['Substitutions%'] = 100 * final_df['Substitutions'] / final_df['Reads_aligned']
    
    # Merge the accession/gene info with amplicon book
    amplicon_book_simple = amplicon_book_long.rename(columns={"fastq": 'Fastq', "PCR_Product_fa": "Amplicon"})
    amplicon_book_simple = amplicon_book_simple[['Fastq', 'Amplicon', 'accession/gene']]
    # Performing the left join
    final_df = pd.merge(final_df, amplicon_book_simple, on=["Fastq", "Amplicon"], how='left')

    # Define the desired column order
    desired_cols_order = [
        'Fastq', 'accession/gene', 'Amplicon', 'Unmodified%', 'Modified%', 
        'Insertions%', 'Deletions%', 'Substitutions%'
    ]

    # Get the remaining columns not specified in the desired order
    remaining_cols = [col for col in final_df.columns if col not in desired_cols_order]

    # Combine the two lists to get the full rearranged order
    cols_order = desired_cols_order + remaining_cols

    # Apply the new column order to the DataFrame
    final_df = final_df[cols_order]
    final_df.to_csv(os.path.join(output_dir, 'summary_editing_frequency.csv'), index=False)
    log_message("Summary of editing frequency saved to 'summary_editing_frequency.csv'")
else:
    log_message("No result files were read successfully.", logging.WARNING)

error_pattern = '**/*error_log.txt'

# Full pattern path combining base directory and the file pattern.
full_error_pattern = os.path.join(output_dir, error_pattern)

# Perform the recursive search to get a list of all matching file paths.
error_file_paths = glob.glob(full_error_pattern, recursive=True)

# Initialize an empty list to store the row data for our future DataFrame.
data = []

# Function to find the line starting with 'ERROR'.
def find_error_line(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('ERROR'):
                return line.strip()  # Return the line without leading/trailing white spaces.
    return None  # Return None if no such line exists.

# Iterate over all the file paths.
for file_path in error_file_paths:
    try:
        # Call the function to get the error line.
        error_line = find_error_line(file_path)

        # Proceed only if an 'ERROR' line was found.
        if error_line:
            # Extract the filename from the full file path.
            base_name = os.path.basename(file_path)
            # Remove the specific part of the file name.
            error_name = base_name.replace('_error_log.txt', '')

            # Append the data to our rows list.
            data.append([error_name, error_line])
    except Exception as e:
        log_message(f"Failed to process file {file_path} due to {e}", logging.ERROR)

# Create a DataFrame from the accumulated data.
if data == []:
    log_message('No errors found.')
else:
    error_df = pd.DataFrame(data, columns=['File_Name', 'Error_Line'])
    error_df.to_csv(os.path.join(output_dir, 'summary_error_log.csv'), index=False)
    log_message("Summary of errors saved to 'summary_error_log.csv'")

log_message("CRISPResso analysis completed")