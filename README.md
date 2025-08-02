# RADseq *Cryophytum crystallinum*

## 🧬 **Linux and R scripts for a RADseq analysis pipeline**

Contact Clarke van Steenderen at vsteenderen@gmail.com or clarke.vansteenderen@ru.ac.za for queries
<br><br> 

<img src="crystallinum.png" width="500">

*Image credit: Centre for Biological Control (CBC)*

## SETUP

The folder structure of each individual project should resemble this, adapted to the number of plates present:

```plaintext
RADseq_project/
└── job_files/
      	├── index_refgenome.job	
      	├── subsample.job
      	├── fastqc.job		
      	├── 1_demultiplex.job
      	├── 2_align_and_stacks_denovo.job	
      	├── 2_align_and_stacks_refgenome.job		
└── data/  
    ├── plate_1/  
    │   ├── filenameA_R1_001.fastq.gz  
    │   └── filenameA_R2_001.fastq.gz  
    └── plate_2/  
    │   ├── filenameB_R1_002.fastq.gz  
    │   └── filenameB_R2_002.fastq.gz  
    └── barcodes/  
    │   ├── internal_indexes_plate_1.txt  
    │   ├── internal_indexes_plate_2.txt  
    │   └── pops_all.txt
    └──ref_genome/	
	├── ncbi_dataset/ (other applicable folders here, if a reference genome is available)
```

On an HPC, create the necessary folders (example below for 2 plates):

```
mkdir -p RADseq_project/data/plate_1 RADseq_project/data/plate_2 RADseq_project/data/barcodes RADseq_project/data/refgenome RADseq_project/jobfiles
```

* Move the job files from this repo into the **jobfiles/** folder
* Move your data from the sequencing company (Read 1 and Read 2 files) into each corresponding plate folder
* If a reference genome is available, download it from NCBI and move it to the **ref_genome/** folder
* Create a sample sheet for the samples on each plate, with the columns below (well, origin, and sample_id are essential).

For example, **samples_plate_1.csv**:

| well | origin       | sample_id | lat      | lon      |
|------|--------------|-----------|----------|----------|
| C2   | Dwarskersbos | DWSA1     | -32.7032 |  18.2203 |
| D3   | Dwarskersbos | DWSA3     | -32.7032 |  18.2203 |
| G6   | Dwarskersbos | DWSA4     | -32.7032 |  18.2203 |
| E9   | Paternoster  | PNSA1     | -32.8096 | 17.89104 |

**samples_plate_2.csv**:

| well | origin       | sample_id | lat      | lon      |
|------|--------------|-----------|----------|----------|
| C1   | Dwarskersbos | DWSA2     | -32.7032 | 18.2203  |
| F3   | Cyprus       | CYP1      | 35.1264  | 33.4299  |
| H3   | Cyprus       | CYP4      | 35.1264  | 33.4299  |
| A4   | Mexico       | MX01.2    | 27.49364 | -114.149 |
| D4   | Mexico       | MX01.3    | 27.45517 | -114.136 |

* If there is more than one plate, create another single file with all sample info combined -> e.g. **all_plates.csv**

| well | origin       | sample_id | lat      | lon      |
|------|--------------|-----------|----------|----------|
| C2   | Dwarskersbos | DWSA1     | -32.7032 |  18.2203 |
| D3   | Dwarskersbos | DWSA3     | -32.7032 |  18.2203 |
| G6   | Dwarskersbos | DWSA4     | -32.7032 |  18.2203 |
| E9   | Paternoster  | PNSA1     | -32.8096 | 17.89104 |
| C1   | Dwarskersbos | DWSA2     | -32.7032 | 18.2203  |
| F3   | Cyprus       | CYP1      | 35.1264  | 33.4299  |
| H3   | Cyprus       | CYP4      | 35.1264  | 33.4299  |
| A4   | Mexico       | MX01.2    | 27.49364 | -114.149 |
| D4   | Mexico       | MX01.3    | 27.45517 | -114.136 |

* Create an internal index file for each plate, and a single population info file for all samples. Use the functions in the R script in the **generate_barcodes/** folder:

```{r}
# source the functions
source("generate_barcodes/barcode_functions.R")

# read in sample info for plate 1
plate_1 = read.csv("samples_plate_1.csv")
head(plate_1)

plate_2 = read.csv("ssamples_plate_2.csv")
head(plate_2)

# run function
plate_1_indexes = create_internal_indexes(plate_samples = plate_1)
plate_2_indexes = create_internal_indexes(plate_samples = plate_2)

# write out as text files
write.table(plate_1_indexes,
           file = "internal_indexes_plate_1.txt", sep = "\t",             
           row.names = FALSE, col.names = FALSE, quote = FALSE, eol = "\n")

write.table(plate_2_indexes,
            file = "internal_indexes_plate_2.txt", sep = "\t",             
            row.names = FALSE, col.names = FALSE, quote = FALSE, eol = "\n")

```

Now create a full population info file for samples across all plates:

```{r}
sample_sheet = read.csv("all_plates.csv")
pops = create_pop_file(sample_sheet = sample_sheet)

# save as a txt file
write.table(pops, file = "sample_info/pops_all.txt", sep = "\t",             
            row.names = FALSE, col.names = FALSE, quote = FALSE, eol = "\n")
```

* Move the internal_indexes_plate_n.txt and pops_all.txt files into the **barcodes/** folder

## RUN LINUX JOB SCRIPTS

For each job, just modify the #PBS paramaters in the files to match your HPC platform, resource and walltime requirements, and logfile paths.

### 🔵 index_refgenome.job

If you have a reference genome, run this script to index it (indexing makes it easier to work with further downstream). 

* REFERENCE_INDEX = file path to the reference genome, which should contain a **.fna** file

To run the script:

```
qsub -v REFERENCE_INDEX="/mnt/lustre/users/cvansteenderen/RADseq_crystallinum/data/ref_genome/ncbi_dataset/data/GCA_030267885.1" 0_index_refgenome.job
```

### 🔵 fastqc.job

This runs a quality check on the data received from the sequencer:

* BASE_DIR = the directory housing your data, containing a separate folder per plate. In the directory template example here -> **data/**		
* NUM_PLATES = the number of plates you have

```
qsub -v BASE_DIR="/mnt/lustre/users/cvansteenderen/RADseq_crystallinum/data/Dean/IcePlant.RawData",NUM_PLATES=2 fastqc.job
```

### 🔵 1_demultiplex.job

This step takes in the internal index information for each sample on each plate (the files you created that are now in **barcodes/**), and searches through the sequence data to find all the fragments that belong to each unique sample. It creates an output folder called **stacksoutput/**, and a folder for each plate. E.g. **stacksoutput/plate_1** and **stacksoutput/plate_2**. Each plate folder will contain a folder (.fq.gz) per sample. Once it has completed demultiplexing, the script creates a new folder called **combined_plates/**, into which it puts all samples from all plates. It then checks for sample folders that are abnormally small, and moves those into a new folder called **removed_zipped/**. The remaining good samples are put into a folder called **ready/**. The barcodes folder is updated to include these ready samples, as the removed ones should no longer be in the sample list. This sample list is saved as **barcodes/bothplates_pops.txt**.		
The file structure should resemble: **data/stacksoutput/combined_plates/ready** and **data/barcodes/bothplates_pops.txt**		
The demultiplexed samples in the **ready/** folder will be paired end reads, denoted by **.1** and **.2** after each sample name. For example, samples CYP1 and DWSA1 would appear as:		

```
CYP1.1.fq.gz				
CYP1.2.fq.gz					
DWSA1.1.fq.gz				
DWSA1.2.fq.gz
```			

Tweak the script to change the Stacks parameters (e.g. enzymes).

Script parameters:

* BASE_DIR = the directory housing your data, containing a separate folder per plate. In the directory template example here -> **data/**		
* NUM_PLATES = the number of plates you have

To run the script:

```
qsub -v NUM_PLATES=2, BASE_DIR="/mnt/lustre/users/cvansteenderen/RADseq_crystallinum/data/Dean_all_spp/IcePlant.RawData" 1_demultiplex.job
```

### 🔵 subsample.job

If your demultiplexed sample files are very large (500MB - 1GB), consider subsampling them (selecting only a fraction of the fragments) before continuing with Stacks. This uses the reformat.sh function in the bbmap library, where N fragments are randomly selected, while keeping the correct Read1 and Read2 pairs.

Script parameters: 

* INPUT_DIR = the path to your ready samples that have been demultiplexed				
* SAMPLERATE = the proportion of fragments to keep	

To run the script:

```
qsub INPUT_DIR="/mnt/lustre/users/cvansteenderen/RADseq_crystallinum/data/Admera/stacksoutput/combined_plates/ready",SAMPLERATE=0.1 subsample.job
```

### ⚠️ Check ⚠️: Are you working with data from more than one project (independent sequence results received from the same provider or different providers at different times)?	
If so: at this point you need to create a new folder in your project's home directory (e.g. **combined_data/ready**), and copy over the sample files from each **ready/** folder into that. In the example below, we have two independent projects that were demultiplexed. Each project path (**ready/** or **ready/subsampled** folder containing final demultiplexed sample files) and its **barcodes/** folder should be listed in the folders and barcodes arrays below. Ignore this if you are only working with data for one project.

```
# Define your HOME_DIR (the base directory of your project)
HOME_DIR="/mnt/lustre/users/cvansteenderen/RADseq_crystallinum"

# Create output directories for the new folder that will contain all samples across projects, and the associated barcodes folder
mkdir -p "$HOME_DIR/combined_data/ready"
mkdir -p "$HOME_DIR/combined_data/barcodes"

# Define arrays of folder and barcode paths (same order). Here, we have two separate project folders with demultiplexed samples in each
# List all the project folders here, pointing to the folder containing the samples that are ready to be processed further
folders=(
  "$HOME_DIR/data/Admera/stacksoutput/combined_plates/ready/subsampled"
  "$HOME_DIR/data/Dean/IcePlant.RawData/stacksoutput/combined_plates/ready"
)

# in the same order, list the associated barcodes folders pointing to the bothplates_pops.txt files -> these contain all the sample info (sample id and pop assignment) for that project
barcodes=(
  "$HOME_DIR/data/Admera/barcodes/bothplates_pops.txt"
  "$HOME_DIR/data/Dean/IcePlant.RawData/barcodes/bothplates_pops.txt"
)

# Loop through folders and copy files over to the new combined_data folder
for folder in "${folders[@]}"; do
  echo "Copying files from: $folder"
  cp "$folder"/* "$HOME_DIR/combined_data/ready/"
done

# Loop through barcode files and concatenate
output_barcode_file="$HOME_DIR/combined_data/barcodes/bothplates_pops.txt"
> "$output_barcode_file"  # Empty or create new

for barcode_file in "${barcodes[@]}"; do
  echo "Appending barcode file: $barcode_file"
  cat "$barcode_file" >> "$output_barcode_file"
done

echo "Combination complete. Output in: $HOME_DIR/combined_data"

```

### 🔵 2_align_and_stacks_denovo.job FOR DENOVO assembly (no reference genome)

This script assembles the demultiplexed paired-end samples using Stacks's **denovo_map.pl** and **populations** workflows. 
It creates an output folder called **stacksoutput_denovo/** and **stacksoutput_denovo/populations/** inside the SAMPLE_DIR. The **populations/** folder will contain a **populations.snps.vcf file** -> this is what you need for downstream SNP analyses.

Check the script to modify Stacks parameters.

Script parameters:

* SAMPLE_DIR = the path to the folder containing the demultiplexed samples that are ready to go (ready/ or ready/subsampled)					
* BARCODES_DIR = the path to the barcodes folder containing the bothplates_pops.txt file

To run the script:

```
qsub -v SAMPLE_DIR="/mnt/lustre/users/cvansteenderen/RADseq_crystallinum/data/combined_data_all_spp/ready",BARCODES_DIR="/mnt/lustre/users/cvansteenderen/RADseq_crystallinum/data/combined_data_all_spp/barcodes" 2_align_and_stacks_denovo.job
```

### 🔵 2_align_and_stacks_refgenome.job FOR ASSEMBLY USING A REFERENCE GENOME

This script differs from the denovo approach in that it aligns sample reads to a reference genome, sorts them, and then assembles the reads using Stacks' **ref_map.pl** function. This script creates an output folder called **refgenome_alignments/** and **refgenome_alignments/populations/** inside the SAMPLE_DIR. The **populations/** folder will contain a **populations.snps.vcf file** -> this is what you need for downstream SNP analyses.

Script parameters:

* SAMPLE_DIR = the path to the folder containing the demultiplexed samples that are ready to go (ready/ or ready/subsampled)					
* BARCODES_DIR = the path to the barcodes folder containing the bothplates_pops.txt file				
* REFERENCE_INDEX = the path to the indexed reference genome

To run the script:

```
qsub -v SAMPLE_DIR="/mnt/lustre/users/cvansteenderen/RADseq_crystallinum/data/combined_data_all_spp/ready",BARCODES_DIR="/mnt/lustre/users/cvansteenderen/RADseq_crystallinum/data/combined_data_all_spp/barcodes",REFERENCE_INDEX="/mnt/lustre/users/cvansteenderen/RADseq_crystallinum/data/ref_genome/ncbi_dataset/data/GCA_030267885.1/reference_index" 2_align_and_stacks_refgenome.job
```

