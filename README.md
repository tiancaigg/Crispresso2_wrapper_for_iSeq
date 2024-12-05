# Crispresso2_wrapper_for_iSeq

A Python wrapper script for running CRISPResso analysis on iSeq data.

## Dependencies

The following packages are required:
- Python 3.x
- pandas
- BioPython
- CRISPResso2
- subprocess
- argparse
- glob
- multiprocessing

### Recommended Installation

We recommend using a conda virtual environment:

```bash
conda create -n crispresso_env python=3.8
conda activate crispresso_env
conda install -c bioconda crispresso2
conda install pandas biopython

##Usage
./Crispresso2_wrapper_v4.py \
    --fastq_dir <path_to_fastq_files> \
    --amplicon_fasta_dir <path_to_fasta_files> \
    --amplicon_book <path_to_excel_booking_list> \
    --sheet <sheet_name> \
    --output_dir <output_directory>


##Example
./Crispresso2_wrapper_v4.py --fastq_dir /media/HDD2/nas/iSeq_fastq_data/20240925_FS10001210_27_BWB90206-1517/Alignment_1/20240926_125353/Fastq/ --amplicon_fasta_dir /home/jinge/iseq_analysis/PCR_fasta/ --amplicon_book /home/jinge/iseq_analysis/iSeq_Amplicon_booking_list_2024_Oct_14.xlsx --sheet Oct2024_run --output_dir test4

##Arguments

--fastq_dir: Directory containing FASTQ files which is already demultiplexed (required).
--amplicon_fasta_dir: Directory containing amplicon FASTA files (with extention ".fa"). name should match the excel file.  (required)
--output_dir: Output directory for CRISPResso (default: ./Crispresso_output)
--amplicon_book: Amplicon booking list in xlsx format (required). includes the registration of amplicon samples and CRISPR gRNA seq.
--sheet: Sheet name in the Amplicon booking list (required)
--verbose: Verbose output (T/F, default: T)

##Output
The script generates several output files:

CRISPResso analysis results for each sample
summary_editing_frequency.csv: Summary of editing frequencies
summary_error_log.csv: Summary of any errors encountered
crispresso_wrapper.log: Detailed log file


