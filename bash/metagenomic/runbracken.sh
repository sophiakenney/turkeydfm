#!/bin/bash
#SBATCH --job-name=runbracken
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=250GB
#SBATCH -t 48:00:00
#SBATCH -o runbracken.out
#SBATCH -e runbracken.err
#SBATCH --export=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user smk459@psu.edu

#--- write job start time to log file
echo "Job started at $(date)"

#--- activate bioinfo env for numpy
source activate base
module load anaconda3
conda activate bioinfo


#--- Set the directory containing the .fastq files
DIR="/storage/group/evk5387/default/sophia_kenney/turkey/shotgun/resistome_v2/reads_decontam/" # set to directory of fasta files
cd ${DIR}

#--- Set variable paths
KRAKEN_DB="/storage/group/evk5387/default/sophia_kenney/kraken_bracken/kraken2/krakendbstandard"
OUTPUT="/storage/home/smk459/scratch/turkey-sg/taxonomy" # set to directory where youd like your output
KRAKEN2_DIR="/storage/group/evk5387/default/sophia_kenney/kraken_bracken/kraken2"
BRAKEN_DIR="/storage/group/evk5387/default/sophia_kenney/kraken_bracken/Bracken"
KRAKEN_TOOLS="/storage/group/evk5387/default/sophia_kenney/kraken_bracken/KrakenTools"
KRONA_TOOLS="/storage/group/evk5387/default/sophia_kenney/kraken_bracken/Krona/KronaTools"


#--- Loop over all .fastq files
for file in "$DIR"*R1.fastq.gz; do

    fname=$(basename $file)
    filename="${fname%.R*}"

    echo "${filename} running"
    echo " "


   # Start Bracken Analysis
    echo "Start Bracken Analysis"
    ${BRAKEN_DIR}/bracken -d ${KRAKEN_DB} -i "${filename}.k2report" -o "${filename}.bracken" -r 151 -w "${filename}.breport" &
    wait $!
    echo "Finished Bracken Analysis"
    echo " "

    # Replace the paths with the path to each specific metric script

    echo "compute alpha diversity"
    python ${KRAKEN_TOOLS}/DiversityTools/alpha_diversity.py -f "${filename}.bracken" -a BP > "${filename}_BP.txt" &
    wait $!
    python ${KRAKEN_TOOLS}/DiversityTools/alpha_diversity.py -f "${filename}.bracken" -a Sh > "${filename}_Sh.txt" &
    wait $!
    python ${KRAKEN_TOOLS}/DiversityTools/alpha_diversity.py -f "${filename}.bracken" -a F > "${filename}_F.txt" &
    wait $!
    python ${KRAKEN_TOOLS}/DiversityTools/alpha_diversity.py -f "${filename}.bracken" -a Si > "${filename}_Si.txt" &
    wait $!
    python ${KRAKEN_TOOLS}/DiversityTools/alpha_diversity.py -f "${filename}.bracken" -a ISi > "${filename}_ISi.txt" &
    wait $!
    echo "alpha diversity complete"

    echo "compute kreport2krona"
    python ${KRAKEN_TOOLS}/kreport2krona.py -r "${filename}.breport" -o "${filename}.b.krona.txt" --no-intermediate-ranks &
    wait $!
    perl ${KRONA_TOOLS}/scripts/ImportText.pl "${filename}.b.krona.txt" -o "${filename}.krona.html" &
    wait $!

    echo "${filename} complete"
    echo " "
done

echo "combine kreports"
python ${KRAKEN_TOOLS}/combine_kreports.py -r ${DIR}*.breport -o "combinedkreports.breport" --display-headers &
wait $!
echo "complete kreport all"
echo " "

echo "compute kreport2krona"
python ${KRAKEN_TOOLS}/kreport2krona.py -r "combinedkreports.breport" -o "combinedkreports.b.krona.txt" &
wait $!
perl ${KRONA_TOOLS}/scripts/ImportText.pl "combinedkreports.b.krona.txt" -o "combinedkreports.krona.html" &
wait $!
echo "complete kreport2krona all"
echo " "

echo "compute beta diversity"
python ${KRAKEN_TOOLS}/DiversityTools/beta_diversity.py -i ${DIR}*.bracken --type bracken > "betadiversity.txt"
echo "completed beta diversity"


# --- move to output directory
echo "start move files to output directory"

mv ${DIR}*.kraken2 ${OUTPUT}
mv ${DIR}*.b.krona.txt ${OUTPUT}
mv ${DIR}*.breport ${OUTPUT}
mv ${DIR}*.k2report ${OUTPUT}
mv ${DIR}*.bracken ${OUTPUT}
mv ${DIR}*.krona.html ${OUTPUT}
mv ${DIR}*.txt ${OUTPUT}
mv ${DIR}*.classified-out.R* ${OUTPUT}/classified
mv ${DIR}*.unclassified-out.R* ${OUTPUT}/unclassified

echo "completed move files to output directory"

#compress fastq files

echo "start compressing classified/unclassified fastq files"

cd ${OUTPUT}

gzip -v ./*classified/*.fastq

echo "completed compressing classified/unclassified fastq files"

#--- write job end time to log file
echo "Job ended at $(date)"
