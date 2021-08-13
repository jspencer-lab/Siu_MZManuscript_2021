#!/usr/bin/env bash
#Clean up and align dataset to create final repertoire fastq file
#Reference https://presto.readthedocs.io/en/stable/workflows/VanderHeiden2017_Workflow.html
##CONNECTING TO DOCKER
#cd C:\Users\Immcantation pipeline\data
#docker run -it -v "C:\Users\Immcantation pipeline\data:/data":z immcantation/suite:4.2.0 bash
for filename in ./*.fastq; do 
    echo "mainloop" ${filename}
    outname=${filename%.fastq}
    outname=${outname:2}
    echo $outname
#From trimmomatic output (removing reads shorter than 150, lower than q 20, remove TruSeq3 adapter sequences), convert to proper presto headings. QC with FastQC.
ConvertHeaders.py illumina -s results/Trimmomatic_paired_V2/output_${outname}_paired.fastq --outdir results/prestoV2/convertheaders
#Remove index
ParseHeaders.py delete -s results/prestoV2/convertheaders/output_${outname}_paired_convert-pass.fastq -f INDEX --outdir results/prestoV2/convertheaders
done

##### NEED TO CHANGE OUTNAME #####
for filename in ./results/prestoV2/assemblepairs/*.fastq; do  
    echo "mainloop" ${filename}
    outname=${filename%_assemble-pass.fastq}
    outname=${outname:33}
    echo $outname
#Assemble R1 and R2 (RC at the tail). Fill gaps using reference (downloaded from Kleinstein)
AssemblePairs.py sequential -1 results/prestoV2/convertheaders/output_${outname}_R2_paired_convert-pass_reheader.fastq \
    -2 results/prestoV2/convertheaders/output_${outname}_R1_paired_convert-pass_reheader.fastq -r IMGT_Human_IG_V.fasta \
    --coord presto --rc tail \
    --aligner blastn --outdir results/prestoV2/assemblepairs --outname ${outname} --log results/prestoV2/assemblepairs/${outname}_AP.log    

#Final processing before IMGT (filter for length, determine C region, collapse duplicates, convert to fasta)
FilterSeq.py length -s results/prestoV2/assemblepairs/${outname}_assemble-pass.fastq -n 200 --outdir results/prestoV2/filterseq --outname ${outname} --log results/prestoV2/filterseq/${outname}_FS.log
ParseLog.py -l results/prestoV2/filterseq/${outname}_FS.log -f ID QUALITY

MaskPrimers.py align -s results/prestoV2/filterseq/${outname}_length-pass.fastq \
     -p AbSeq_Human_IG_InternalCRegion.fasta --maxlen 100 --maxerror 0.3 \
     --mode tag --revpr --skiprc --pf CREGION --outdir results/prestoV2/maskprimers2 --outname ${outname} --log results/prestoV2/maskprimers2/${outname}_MP3.log 

#Collapse duplicate sequences (must have same C region). also removes sequences with a high number of interior N-valued nucleotides (-n 20 and --inner)
CollapseSeq.py -s results/prestoV2/maskprimers2/${outname}_primers-pass.fastq -n 20 --inner \
     --uf CREGION --outname ${outname} --outdir results/prestoV2/collapseseq --log results/prestoV2/collapseseq/${outname}_collapseseq.log 

seqtk seq -a results/prestoV2/collapseseq/${outname}_collapse-unique.fastq > results/prestoV2/collapseseq_fasta/${outname}_collapse-unique.fasta
done 

#After IMGT-heavy run (default setting, AIRR format)
for filename in ./results/prestoV2/IMGT/*.tsv; do 
    echo "mainloop" ${filename}
    outname=${filename%_unique_db-pass.tsv}
    outname=${outname:9}
    echo $outname

    #Parsing IMGT results
    MakeDb.py imgt -i results/prestoV2/IMGT/22072021_${outname}_unique.txz -s results/prestoV2/collapseseq_fasta/${outname}_collapse-unique.fasta --extended

    #Filtering for productive sequences only
    ParseDb.py select -d ${filename} -f productive -u T
done

#After determining clonal threshold in Shazam, merge tissues, adding clone definitions (same V-J, same length, similar distance) in R
DefineClones.py -d results/prestoV2/IMGT/donorA-alltissues_parse-select.tsv --act set --model ham --norm len --dist 0.1 --format airr --nproc 16 --outname donorA_combined --outdir results/prestoV2/IMGT/
DefineClones.py -d results/prestoV2/IMGT/donorB-alltissues_parse-select.tsv --act set --model ham --norm len --dist 0.11 --format airr --nproc 16 --outname donorB_combined --outdir results/prestoV2/IMGT/
DefineClones.py -d results/prestoV2/IMGT/donorC-alltissues_parse-select.tsv --act set --model ham --norm len --dist 0.13 --format airr --nproc 16 --outname donorC_combined --outdir results/prestoV2/IMGT/

#After combining with bulk and sc data in R, re-define clones (same V-J, same length, similar distance)
DefineClones.py -d results/combo/donorA-combo_dat.tsv --act set --model ham --norm len --dist 0.09 --format airr --nproc 16 --outname donorA_datcombined09 --outdir results/combo
DefineClones.py -d results/combo/donorB-combo_dat.tsv --act set --model ham --norm len --dist 0.11 --format airr --nproc 16 --outname donorB_datcombined11 --outdir results/combo
DefineClones.py -d results/combo/donorC-combo_dat.tsv --act set --model ham --norm len --dist 0.13 --format airr --nproc 16 --outname donorC_datcombined --outdir results/combo
