#!/bin/bash
#SBATCH --job-name=runkraken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=250GB
#SBATCH -t 48:00:00
#SBATCH -o runkraken.out
#SBATCH -e runkraken.err
#SBATCH --export=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user smk459@psu.edu


# Kraken 2 has two output files, a standard Kraken 2 output,
# .kraken2 file, and a Kraken 2 report,.k2report file. The
# standard Kraken 2 file outputs lines containing five tab-delimited fields;
# (1) ‘C’/‘U’ indicating classified or unclassified,
# (2) sequence ID,
# (3) taxonomy ID (0 if unclassified),
# (4) length of sequence in base pairs,
# (5) a space-delimited list indicating the lowest common ancestor mapping of each k-mer.

# The Kraken report file outputs lines containing six tab-delimited fields:
# (1) percentage of fragments covered by the clade rooted at this taxon,
# (2) number of fragments covered by the clade rooted at this taxon,
# (3) number of fragments assigned directly to this taxon,
# (4) rank code such as S for species,
# (5) NCBI taxonomic ID,
# (6) indented scientific name.

#--- write job start time to log file
echo "Job started at $(date)"

#--- Set the directory containing the .fastq files
DIR="/storage/group/evk5387/default/sophia_kenney/turkey/shotgun/resistome_v2/reads_decontam/" # set to directory of fasta files
cd ${DIR}

#--- Set variable paths
KRAKEN_DB="/storage/group/evk5387/default/sophia_kenney/kraken_bracken/kraken2/krakendbstandard"
OUTPUT="/storage/home/smk459/scratch/kraken-out/" # set to directory where youd like your output
KRAKEN2_DIR="/storage/group/evk5387/default/sophia_kenney/kraken_bracken/kraken2"

#--- Loop over all .fq files
for file in "$DIR"*R1.fastq.gz; do

    fname=$(basename $file)
    filename="${fname%.R*}"

    echo "${filename} running"
    echo " "

    echo "kraken2 2"
    ${KRAKEN2_DIR}/kraken2 --threads 16 --db ${KRAKEN_DB} \
    --gzip-compressed \
    --paired --classified-out ${filename}.classified-out.R#.fastq \
    --unclassified-out ${filename}.unclassified-out.R#.fastq \
    --report "${filename}.k2report" \
    --report-minimizer-data \
    --minimum-hit-groups 3 \
    ${filename}.R1.fastq ${filename}.R2.fastq > "${filename}.kraken2" &
    wait $!
    echo "kraken2 2 complete"
done

#--- write job end time to log file
echo "Job ended at $(date)"
