# Cladocopium goreaui Transcriptome Annotation, version December 18, 2024
# Created by Misha Matz (matz@utexas.edu), modified by Michael Studivan (studivanms@gmail.com) for use on FAU's HPC (KoKo)


#------------------------------
# BEFORE STARTING, replace, in this whole file:
#	- studivanms@gmail.com by your actual email;
#	- mstudiva with your KoKo user name.

# The idea is to copy the chunks separated by empty lines below and paste them into your cluster
# terminal window consecutively.

# The lines beginning with hash marks (#) are explanations and additional instructions -
# please make sure to read them before copy-pasting.


#------------------------------
# setup

# To install Bioperl as a conda environment
conda create -y -n bioperl perl-bioperl

# getting scripts
cd ~/bin
git clone https://github.com/z0on/annotatingTranscriptomes.git
mv annotatingTranscriptomes/* .
rm -rf annotatingTranscriptomes
rm launcher_creator.py

git clone https://github.com/z0on/emapper_to_GOMWU_KOGMWU.git
mv emapper_to_GOMWU_KOGMWU/* .
rm -rf emapper_to_GOMWU_KOGMWU

git clone https://github.com/mstudiva/Cladocopium-annotated-transcriptome.git
mv Cladocopium-annotated-transcriptome/* .
rm -rf Cladocopium-annotated-transcriptome

# creating backup directory
mkdir backup

# creating annotation directory
cd
mkdir annotate
cd annotate


#------------------------------
# getting transcriptomes

# Chen (2022) https://doi.org/10.3390/microorganisms10081662
# cds and protein translations from genome at https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_947184155.2/

# Renaming gene identifiers for ease
sed -i 's/lcl|CAMXCT/Cladocopium/' Cladocopium.fasta
sed -i 's/C1SCF055_LOCUS/Cladocopium/' Cladocopium.fasta

# transcriptome statistics
conda activate bioperl
echo "seq_stats.pl Cladocopium.fasta > seqstats_Cladocopium.txt" > seq_stats
launcher_creator.py -j seq_stats -n seq_stats -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch seq_stats.slurm

Cladocopium.fasta
-------------------------
49545 sequences.
2003 average length.
99816 maximum length.
129 minimum length.
N50 = 2910
99.3 Mb altogether (99250232 bp).
0 ambiguous Mb. (0 bp, 0%)
0 Mb of Ns. (0 bp, 0%)
-------------------------


#-------------------------
# Extracting contig, isogroup, and protein IDs into lookup tables

grep ">" Cladocopium.fasta | awk '{
    header = substr($1, 2);                  # Remove the leading ">" from the first field
    match($0, /locus_tag=([^]]+)/, gene);    # Extract the locus_tag (gene name)
    print header "\t" gene[1];               # Print the full header and locus_tag
}' > Cladocopium_seq2iso.tab

grep ">" Cladocopium.fasta | awk -F'[][]' '{
    for (i=1; i<=NF; i++) {
        if ($i ~ /locus_tag=/) { gsub("locus_tag=", "", $i); gene=$i }
        if ($i ~ /protein_id=/) { gsub("protein_id=", "", $i); protein=$i }
    }
    if (gene && protein) { print protein, gene }
}' > Cladocopium_seq2pro.tab


#------------------------------
# GO annotation
# updated based on Misha Matz's new GO and KOG annotation steps on GitHub: https://github.com/z0on/emapper_to_GOMWU_KOGMWU

# first rename protein headers as gene ID based on lookup table
awk 'BEGIN {
    while (getline < "Cladocopium_seq2pro.tab") {
        map[$1] = $2
    }
}
/^>/ {
    protein_id = substr($0, 2, index($0, " ") - 2)
    if (protein_id in map) {
        sub(protein_id, map[protein_id])
    }
}
{ print }' protein.faa > Cladocopium_pro.fasta

# selecting the longest contig per isogroup
fasta2SBH.pl Cladocopium_pro.fasta >Cladocopium_pro_out.fasta

# scp *_pro_out.fasta to laptop, submit it to
http://eggnog-mapper.embl.de
cd /path/to/local/directory
scp mstudiva@koko-login.hpc.fau.edu:~/path/to/HPC/directory/\*_pro_out.fasta .

# copy link to job ID status and output file, paste it below instead of current link:
# http://eggnog-mapper.embl.de/job_status?jobname=MM_a3tl4oxj

# once it is done, download results to HPC:
wget http://eggnog-mapper.embl.de/MM_a3tl4oxj/out.emapper.annotations

# GO:
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$10 }' out.emapper.annotations | grep GO | perl -pe 's/,/;/g' >Cladocopium_iso2go.tab
# gene names:
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$8 }' out.emapper.annotations | grep -Ev "\tNA" >Cladocopium_iso2geneName.tab


#------------------------------
# KOG annotation
# updated based on Misha Matz's new GO and KOG annotation steps on github: https://github.com/z0on/emapper_to_GOMWU_KOGMWU

cp ~/bin/kog_classes.txt .

#  KOG classes (single-letter):
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$7 }' out.emapper.annotations | grep -Ev "\tNA" >Cladocopium_iso2kogClass1.tab
# converting single-letter KOG classes to text understood by KOGMWU package (must have kog_classes.txt file in the same dir):
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1] = $2;next} {print $1,a[$2]}' kog_classes.txt Cladocopium_iso2kogClass1.tab > Cladocopium_iso2kogClass.tab


#------------------------------
# KEGG annotations:

# first rename fasta headers as gene ID (isogroup) rather than contig ID
awk 'BEGIN {
    while (getline < "Cladocopium_seq2iso.tab") {
        map[$1] = $2
    }
}
/^>/ {
    gene_id = substr($0, 2, index($0, " ") - 2)
    if (gene_id in map) {
        sub(gene_id, map[gene_id])
    }
}
{ print }' Cladocopium.fasta > Cladocopium_iso.fasta

# selecting the longest contig per isogroup
fasta2SBH.pl Cladocopium_iso.fasta >Cladocopium_iso_out.fasta

# Sanity check: How many unique isogroups do we have?
grep ">" Cladocopium.fasta | sort | uniq | wc -l            # 49545
grep ">" Cladocopium_pro.fasta | sort | uniq | wc -l        # 45371
grep ">" Cladocopium_pro_out.fasta | sort | uniq | wc -l    # 45316
grep ">" Cladocopium_iso.fasta | sort | uniq | wc -l        # 49545
grep ">" Cladocopium_iso_out.fasta | sort | uniq | wc -l    # 45319
# The two _out files should roughly match

# scp *_iso_out.fasta to your laptop
cd /path/to/local/directory
scp mstudiva@koko-login.hpc.fau.edu:~/path/to/HPC/directory/\*_iso_out.fasta .
# use web browser to submit _iso.fasta file to KEGG's KAAS server (http://www.genome.jp/kegg/kaas/)
# select SBH method, upload nucleotide query
https://www.genome.jp/kaas-bin/kaas_main?mode=user&id=1734546497&key=banMjCtB

# Once it is done, download to HPC - it is named query.ko by default
wget https://www.genome.jp/tools/kaas/files/dl/1734546497/query.ko

# selecting only the lines with non-missing annotation:
cat query.ko | awk '{if ($2!="") print }' > Cladocopium_iso2kegg.tab

# the KEGG mapping result can be explored for completeness of transcriptome in terms of genes found,
# use 'html' output link from KAAS result page, see how many proteins you have for conserved complexes and pathways, such as ribosome, spliceosome, proteasome etc


#------------------------------
# file transfer

# copy all files to local machine
cd /path/to/local/directory
scp mstudiva@koko-login.fau.edu:~/path/to/HPC/directory/\* .
