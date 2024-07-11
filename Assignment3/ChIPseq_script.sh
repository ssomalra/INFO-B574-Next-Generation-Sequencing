##### Load modules 
module load samtools
module load fastqc
module load trimgalore
module load bowtie/2.5.1
#module load picard
pip install macs3 --user
#pip install deeptools --user # if not installed

mkdir -p /N/scratch/ssomalra/ChIPseq
cd /N/scratch/ssomalra/ChIPseq

##Download ngsutilsj
mkdir -p Apps
wget https://github.com/compgen-io/ngsutilsj/releases/download/ngsutilsj-0.5.0/ngsutilsj -P Apps
chmod +x Apps/ngsutilsj
wget https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar -P Apps
picard=/N/scratch/ssomalra/ChIPseq/Apps/picard.jar


##download demo data
mkdir Data
wget http://tracks.ccbb.iupui.edu/tracks/sliu19/NGS2024/Control_R1_L001.fastq.gz -P Data/
wget http://tracks.ccbb.iupui.edu/tracks/sliu19/NGS2024/Control_R2_L001.fastq.gz -P Data/

wget http://tracks.ccbb.iupui.edu/tracks/sliu19/NGS2024/H3K4me3_S1_R1_L001.fastq.gz -P Data/
wget http://tracks.ccbb.iupui.edu/tracks/sliu19/NGS2024/H3K4me3_S1_R2_L001.fastq.gz -P Data/

## download Index files
wget http://tracks.ccbb.iupui.edu/tracks/sliu19/NGS2024/genome.1.bt2 -P Index
wget http://tracks.ccbb.iupui.edu/tracks/sliu19/NGS2024/genome.2.bt2 -P Index
wget http://tracks.ccbb.iupui.edu/tracks/sliu19/NGS2024/genome.3.bt2 -P Index
wget http://tracks.ccbb.iupui.edu/tracks/sliu19/NGS2024/genome.4.bt2 -P Index
wget http://tracks.ccbb.iupui.edu/tracks/sliu19/NGS2024/genome.rev.1.bt2 -P Index
wget http://tracks.ccbb.iupui.edu/tracks/sliu19/NGS2024/genome.rev.2.bt2 -P Index
wget http://tracks.ccbb.iupui.edu/tracks/sliu19/NGS2024/genome.fa -P Index


SAMPLE="H3K4me3_S1"
CONTROL="Control"

#FASTQC
mkdir -p report/fastqc/
fastqc Data/*.fastq.gz -o report/fastqc/

#adapter trimming and fastqc
mkdir -p report/trimming/${SAMPLE}
trim_galore --paired --fastqc Data/${SAMPLE}_R1_L001.fastq.gz Data/${SAMPLE}_R2_L001.fastq.gz -o report/trimming/${SAMPLE}/
#CONTROL
mkdir -p report/trimming/${CONTROL}
trim_galore --paired --fastqc Data/${CONTROL}_R1_L001.fastq.gz Data/${CONTROL}_R2_L001.fastq.gz -o report/trimming/${CONTROL}/

#alignment and sort
mkdir -p report/alignment/${SAMPLE}
bowtie2 --no-mixed -p 2 -x Index/genome -1 report/trimming/${SAMPLE}/${SAMPLE}_R1_L001_val_1.fq.gz -2 report/trimming/${SAMPLE}/${SAMPLE}_R2_L001_val_2.fq.gz | samtools view -Su /dev/stdin | samtools sort /dev/stdin -o report/alignment/${SAMPLE}/${SAMPLE}_sorted.bam
#CONTROL
mkdir -p report/alignment/${CONTROL}
bowtie2 --no-mixed -p 2 -x Index/genome -1 report/trimming/${CONTROL}/${CONTROL}_R1_L001_val_1.fq.gz -2 report/trimming/${CONTROL}/${CONTROL}_R2_L001_val_2.fq.gz | samtools view -Su /dev/stdin | samtools sort /dev/stdin -o report/alignment/${CONTROL}/${CONTROL}_sorted.bam

#mark duplicates
java -jar $picard MarkDuplicates --INPUT report/alignment/${SAMPLE}/${SAMPLE}_sorted.bam --OUTPUT report/alignment/${SAMPLE}/${SAMPLE}_MD.bam --METRICS_FILE report/alignment/${SAMPLE}/${SAMPLE}_PCRDupes.txt
#CONTROL
java -jar $picard MarkDuplicates --INPUT report/alignment/${CONTROL}/${CONTROL}_sorted.bam --OUTPUT report/alignment/${CONTROL}/${CONTROL}_MD.bam --METRICS_FILE report/alignment/${CONTROL}/${CONTROL}_PCRDupes.txt


#filtering and indexing bam file
Apps/ngsutilsj bam-filter --mapped --no-qcfail --tag-min MAPQ:30 report/alignment/${SAMPLE}/${SAMPLE}_MD.bam report/alignment/${SAMPLE}/${SAMPLE}_MD_filtered.bam
samtools index report/alignment/${SAMPLE}/${SAMPLE}_MD_filtered.bam
#CONTROL
Apps/ngsutilsj bam-filter --mapped --no-qcfail --tag-min MAPQ:30 report/alignment/${CONTROL}/${CONTROL}_MD.bam report/alignment/${CONTROL}/${CONTROL}_MD_filtered.bam
samtools index report/alignment/${CONTROL}/${CONTROL}_MD_filtered.bam


##QC fingerprint
mkdir -p report/QC/${SAMPLE}/
plotFingerprint -b report/alignment/${SAMPLE}/${SAMPLE}_MD_filtered.bam -plot report/QC/${SAMPLE}/${SAMPLE}_fingerprint.png --ignoreDuplicates -l ${SAMPLE}
#CONTROL
mkdir -p report/QC/${CONTROL}/
plotFingerprint -b report/alignment/${CONTROL}/${CONTROL}_MD_filtered.bam -plot report/QC/${CONTROL}/${CONTROL}_fingerprint.png --ignoreDuplicates -l ${CONTROL}


##QC estimate fragment size
bamPEFragmentSize --bamfiles  report/alignment/${SAMPLE}/${SAMPLE}_MD_filtered.bam -hist report/QC/${SAMPLE}/${SAMPLE}_fragment.png --samplesLabel ${SAMPLE}
#CONTROL
bamPEFragmentSize --bamfiles  report/alignment/${CONTROL}/${CONTROL}_MD_filtered.bam -hist report/QC/${CONTROL}/${CONTROL}_fragment.png --samplesLabel ${CONTROL}

#call narrowpeak - narrowPeak, automatically estimate fragmentsize
mkdir -p report/callpeak/${SAMPLE}/
macs3 callpeak -t report/alignment/${SAMPLE}/${SAMPLE}_MD_filtered.bam -c report/alignment/${CONTROL}/${CONTROL}_MD_filtered.bam -f BAM -g hs -q 0.01 -n ${SAMPLE}_MACS3 --outdir report/callpeak/${SAMPLE}

#call broadpeak - broadPeak, automatically estimate fragmentsize
macs3 callpeak -t report/alignment/${SAMPLE}/${SAMPLE}_MD_filtered.bam -c report/alignment/${CONTROL}/${CONTROL}_MD_filtered.bam -f BAM -g hs --broad -n ${SAMPLE}_MACS3_broad --outdir report/callpeak/${SAMPLE}



