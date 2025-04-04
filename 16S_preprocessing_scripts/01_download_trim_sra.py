"""
A script to download (and optionally trim) SRA raw sequence data from an accession list


Inputs ---------------------------------------

DATASET_DIRS
    - A list of directories for each dataset to download
 
ACC_LIST_FILENAME 
    - A csv file which contains accessions for the samples in the dataset
    - One file per dataset will be found in each directory in DATASET_DIRS
    - Usually "SraAccList.csv"
    - This file can be downloaded from the NCBI BioProject site directly 
        - e.g. go to https://www.ncbi.nlm.nih.gov/bioproject/401920
        - click the number of SRA Experiments 
        - and press Send to > File > Format > Accessions List > Create File

        
Requirements ---------------------------------------

SRA Toolkit
    - Tutorial: https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit

Trimmomatic

Cutadapt

"""


# Imports
import os
import subprocess
import csv
import gzip
import shutil



# General Options ==================
# Directories where the datasets are 
DATASET_DIRS = ["D:/microbiome_sequencing_datasets/celiac_16s_datasets/16S_179_Verdu/"]
# The accession list files name (found in these directories ^^^) (.csv or .txt is fine)
ACC_LIST_FILENAME = "SraAccListReduced.csv"
# The subdirectories to put the... 
DOWNLOAD_DIR = "sra_downloads/" # .sra files 
FASTQS_DIR = "fastqs/" # .fastq files
TRIMMED_DIR = "trimmed/" # trimmed .fastq files
# Deletes the SRA, untrimmed, ZIP files when done to save space
DELETE_SRA_AFTER = True # After dumping
DELETE_UNTRIMMED_AFTER = False # After trimming


# Trimming Options =================
# Only trim fwd reads
SINGLE_READS = False
# Delete any outputted unpaired reads after trimming?
DELETE_UNPAIRED_AFTER = True


# Trimmomatic Options =================
# Trim with trimmomatic?
TRIM_WITH_TRIMMOMATIC = False
# Path to trimmomatic
TRIMMOMATIC_PATH = "Trimmomatic-0.39/trimmomatic-0.39.jar"
# Select the trimmomatic adapter file!
# "TruSeq2 is used in GAII machines and TruSeq3 is used by HiSeq and MiSeq machines"
TRIMMOMATIC_ADAPTERS = "Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa"
# Suffixes of dumped fastq files to know which is fw and rv for trimmomatic
FW_READ_SUFFIX = "_1.fastq.gz"
RV_READ_SUFFIX = "_2.fastq.gz"
# Num thrreads
THREADS = 10


# Cutadapt Options =================
# Trim with Cutadapt?
TRIM_WITH_CUTADAPT = True
# Used for 16S_20_Caminero and 16S_179_Verdu:
FWD_CUTADAPT_ADAPTERS = [
    ("GTATTACCGCGGCTGCTGG", 'front')
    ]
REV_CUTADAPT_ADAPTERS = [
    ("CCTACGGGAGGCAGCAG", 'front')
    ]



def trim_reads_trimmomatic_paired(input_forward, input_reverse, accession, output_dir, delete_after=False):

    # Output file paths
    output_forward_paired = os.path.join(output_dir, accession + "_1.fq.gz")
    output_forward_unpaired = os.path.join(output_dir, accession + "_1_unpaired.fq.gz")
    output_reverse_paired = os.path.join(output_dir, accession + "_2.fq.gz")
    output_reverse_unpaired = os.path.join(output_dir, accession + "_2_unpaired.fq.gz")

    # ILLUMINACLIP
    seed_mismatches = "2"
    palindrome_clip_threshold = "30"
    simple_clip_threshold = "10"
    min_adapter_length = "4"
    keep_both_reads = "true"
    illumina_clip_options = f"ILLUMINACLIP:{TRIMMOMATIC_ADAPTERS}:{seed_mismatches}:{palindrome_clip_threshold}:{simple_clip_threshold}:{min_adapter_length}:{keep_both_reads}"

    # Build the command list
    command = [
        "java", "-jar", TRIMMOMATIC_PATH, "PE", "-threads", str(THREADS),
        input_forward, input_reverse,
        output_forward_paired, output_forward_unpaired,
        output_reverse_paired, output_reverse_unpaired,
        illumina_clip_options
    ]

    # Execute the command
    subprocess.call(command)

    # Delete input files
    if delete_after:
        os.remove(input_forward)
        os.remove(input_reverse)

    return output_forward_paired, output_reverse_paired

def unzip_files(file_paths, delete_after=True):
    unzipped_files = []
    for file_path in file_paths:
        if not file_path.endswith('.gz'):
            continue
        output_path = file_path[:-3]  # Remove the .gz extension
        # Unzip the file
        with gzip.open(file_path, 'rb') as f_in:
            with open(output_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        unzipped_files.append(output_path)
        # Delete the original .gz file if requested
        if delete_after:
            os.remove(file_path)
    return unzipped_files

def write_acc_list(accessions, dir=os.getcwd()):
    # Write acc_list.txt file with accession for each line
    with open(os.path.join(dir, "acc_list.txt"), "w") as f:
        for accession in accessions:
            f.write(accession)
            if accession != accessions[-1]:
                f.write("\n")

def read_acc_list(acc_list_file):
    accessions = []
    # If a txt file
    if acc_list_file.endswith(".txt"):
        with open(acc_list_file, "r") as f:
            for line in f:
                accession = line.strip()  # Remove leading/trailing whitespace
                if len(accession) != 0:  # Skip empty lines
                    accessions.append(accession)
    # If a csv file
    elif acc_list_file.endswith(".csv"):
        with open(acc_list_file, 'r', newline='') as csvfile:
            reader = csv.reader(csvfile)
            # Skip the header row
            next(reader)
            for row in reader:
                accessions.append(row[0])
    else:
        raise Exception("Invalid file type")
    return accessions

def download_sras_from_acc_list(download_dir):
    # Download the data using prefetch
    # Run prefetch command
    # os.system('export PATH="/home/haig/sratoolkit.3.1.0-ubuntu64/bin:$PATH"') # Linux
    os.system("prefetch -O " + download_dir + " --option-file acc_list.txt")

def fastq_dump_sras(sra_file_paths, zips_dir, delete_after=False):
    all_zip_paths = []
    # For each sra file path
    for sra_filepath in sra_file_paths:
        # Execute fastq-dump
        subprocess.run(['fastq-dump', '--split-3', '--gzip', "--outdir", zips_dir, sra_filepath])
        if delete_after:
            # Delete the .sra file
            os.remove(sra_filepath)
        # Get the zip file path
        acc = os.path.splitext(os.path.basename(sra_filepath))[0]
        # Get all ".fastq.gz" files in zips_dir with acc in the name
        zip_paths = [zips_dir + filename for filename in os.listdir(zips_dir) if filename.endswith(".fastq.gz") and filename.startswith(acc)]
        # Append to list
        all_zip_paths += zip_paths
    # Return all resulting zip files
    return all_zip_paths

def make_cutadapt_command_single(adapters, in_file):
    # Initialize the base command
    command = "cutadapt"
    
    # Add adapters to the command
    for adapter, position in adapters:
        if position == 'front':
            command += f' -g "{adapter};e=0.25"'
        elif position == 'back':
            command += f' -a "{adapter};e=0.25"'
    
    # Generate the output file name
    out_file = in_file.replace(".fastq.gz", "_trimmed.fastq.gz").replace(".fq.gz", "_trimmed.fq.gz")
    
    # Add the output file, number of rounds, and input file
    command += f" -o {out_file} --minimum-length=10 -n=10 {in_file}"
    
    return command, [out_file]

def make_cutadapt_command_paired(fwd_adapters, fwd_in_file, rev_adapters, rev_in_file):
    # Initialize the base command
    command = "cutadapt"
    
    # Add fwd adapters to the command
    for adapter, position in fwd_adapters:
        if position == 'front':
            command += f' -g "{adapter};e=0.25"'
        elif position == 'back':
            command += f' -a "{adapter};e=0.25"'
    # Add rev adapters to the command
    for adapter, position in rev_adapters:
        if position == 'front':
            command += f' -G "{adapter};e=0.25"'
        elif position == 'back':
            command += f' -A "{adapter};e=0.25"'
    
    # Generate the output file name
    fwd_out_file = fwd_in_file.replace(".fastq.gz", "_trimmed.fastq.gz").replace(".fq.gz", "_trimmed.fq.gz")
    rev_out_file = rev_in_file.replace(".fastq.gz", "_trimmed.fastq.gz").replace(".fq.gz", "_trimmed.fq.gz")
    
    # Add the output file, number of rounds, and input file
    command += f" -o {fwd_out_file} -p {rev_out_file} --minimum-length=10 -n=10 {fwd_in_file} {rev_in_file}"
    
    out_files = [fwd_out_file, rev_out_file]
    return command, out_files

def run_cutadapt_command(command, in_files, delete_after=True):
    try:
        # Execute the command
        subprocess.run(command, shell=True, check=True)
        print(f"Cutadapt command executed successfully for {in_files}")
        
        # Delete input file if specified
        if delete_after:
            for in_file in in_files:
                os.remove(in_file)
                print(f"Deleted input file: {in_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error while executing Cutadapt: {e}")
    except FileNotFoundError as e:
        print(f"File not found: {e}")




# For every dataset
for dataset_dir in DATASET_DIRS:
    # List of accessions
    accessions = read_acc_list(os.path.join(dataset_dir, ACC_LIST_FILENAME))
    # Write it as a txt (if it isnt already)
    if not os.path.exists(os.path.join(dataset_dir, "acc_list.txt")):
        write_acc_list(accessions, dir=dataset_dir)
    print("Will now attempt to download: ", accessions)
    # Ensure this directory exists
    trimmed_dir = os.path.join(dataset_dir, TRIMMED_DIR)
    if not os.path.exists(trimmed_dir):
        os.makedirs(trimmed_dir)
    
    # For each accession
    for i, accession in enumerate(accessions):

        # Download and unpack ==================================================
        print("\n=======================")
        print("Starting process (" + str(i + 1) + "/" + str(len(accessions)) + "): " + accession)

        # Write the single accession to file
        write_acc_list([accession])
        
        # Run prefetch to download SRA files
        print("Attempting to download " + accession + " SRA file...")
        download_sras_from_acc_list(os.path.join(dataset_dir, DOWNLOAD_DIR))
        
        # Look to see if the SRA files downloaded correctly
        sra_paths = [dataset_dir + DOWNLOAD_DIR + accession + "/" + filename for filename in os.listdir(dataset_dir + DOWNLOAD_DIR + accession + "/") if (filename.endswith(".sra") or filename.endswith(".sralite")) and filename.startswith(accession)]
        if len(sra_paths) == 0:
            print("No SRA files detected")
        else:
            print("It appears the SRA files downloaded correctly: ", sra_paths)

        # Run fastq-dump to turn SRAs into zip files (and then delete sras)
        print("\nDumping SRAs: " + str(sra_paths) + "...")
        zip_paths = fastq_dump_sras(sra_paths, dataset_dir + FASTQS_DIR, delete_after=DELETE_SRA_AFTER)
        
        # Look to see if the SRA files dumped correctly
        zip_paths = [dataset_dir + FASTQS_DIR + filename for filename in os.listdir(dataset_dir + FASTQS_DIR) if (filename.endswith(".fq.gz") or filename.endswith(".fastq.gz")) and filename.startswith(accession)]
        if len(zip_paths) == 0:
            print("No dumped SRAs (into zips) detected")
        else:
            print("It appears the SRA files were dumped correctly to zips: ", zip_paths)

        


        # Decide which files are forward and reverse
        input_forward = next((path for path in zip_paths if path.endswith(FW_READ_SUFFIX)), None)
        input_reverse = next((path for path in zip_paths if path.endswith(RV_READ_SUFFIX) and not SINGLE_READS), None)
        if not input_forward:
            raise ValueError(f"Forward reads with suffix {FW_READ_SUFFIX} not detected!!!")
        if not input_reverse and not SINGLE_READS:
            raise ValueError(f"Reverse reads with suffix {RV_READ_SUFFIX} not detected!!! (set SINGLE_READS True if single reads)")




        # Trimmomatic ==================================================
        if TRIM_WITH_TRIMMOMATIC:

            # If single end
            if SINGLE_READS or REV_CUTADAPT_ADAPTERS == []:
                raise ValueError(f"Singles reads not supported.")
        
            # If paired end
            else:
                # Trim paired reads
                output_forward_paired, output_reverse_paired = trim_reads_trimmomatic_paired(input_forward, input_reverse, accession, trimmed_dir, delete_after=DELETE_UNTRIMMED_AFTER)



        # Cutadapt ==================================================
        if TRIM_WITH_CUTADAPT:

            # If single end
            if SINGLE_READS or REV_CUTADAPT_ADAPTERS == []:
                # Cut forward read
                command, out_files = make_cutadapt_command_single(FWD_CUTADAPT_ADAPTERS, input_forward)
                run_cutadapt_command(command, [input_forward], delete_after=DELETE_UNTRIMMED_AFTER)
                
            # If paired end
            else:
                command, out_files = make_cutadapt_command_paired(FWD_CUTADAPT_ADAPTERS, input_forward, REV_CUTADAPT_ADAPTERS, input_reverse)
                run_cutadapt_command(command, [input_forward, input_reverse], delete_after=DELETE_UNTRIMMED_AFTER)


    # If deleting SRA files
    if DELETE_SRA_AFTER:
        shutil.rmtree(dataset_dir + DOWNLOAD_DIR)


# Delete the temporary acc_list.txt
os.remove("acc_list.txt")

# Done
print("Done!")
