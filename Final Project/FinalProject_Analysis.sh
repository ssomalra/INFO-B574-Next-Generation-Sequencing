### Sample Information
# Treatment: ENCFF522PUA and ENCFF094LXX
# Control: ENCFF247MLF and ENCFF086CFB

### Download Data
mkdir fastq_files
wget https://www.encodeproject.org/files/ENCFF522PUA/@@download/ENCFF522PUA.fastq.gz -P fastq_files/treatment1.fastq.gz
wget https://www.encodeproject.org/files/ENCFF094LXX/@@download/ENCFF094LXX.fastq.gz -P fastq_files/treatment2.fastq.gz
wget https://www.encodeproject.org/files/ENCFF247MLF/@@download/ENCFF247MLF.fastq.gz -P fastq_files/control1.fastq.gz
wget https://www.encodeproject.org/files/ENCFF086CFB/@@download/ENCFF086CFB.fastq.gz -P fastq_files/control2.fastq.gz

gunzip fastq_files/*

### Quality Control
module load fastqc
mkdir fastqc_report
fastqc fastq_files/*.fastq -o fastqc_report/

### Trimming
pip install cutadapt
module load trimgalore
trim_galore --phred33 --fastqc fastq_files/*.fastq -o fastq_files/
trim_galore --phred33 --fastqc fastq_files/*.fastq -o fastq_files/

### Alignment
# Get reference genome
wget ftp://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Build index
module load hisat
mkdir assembly
hisat2-build Homo_sapiens.GRCh38.dna.primary_assembly.fa assembly/genome_105

# Alignment
mkdir alignment
mkdir alignment/alignment_report
hisat2 -q -x assembly/genome_105 -U fastq/*_trimmed.fq -S alignment/*.sam --summary-file alignment/alignment_report/*_summary.txt

### Binary Index SAM Files
module load samtools
samtools view -S -b alignment/*.sam > alignment/*.bam
samtools sort alignment/*.bam -o alignment/*.sorted.bam
samtools index -b alignment/*.sorted.bam

### Peak Calling
# NOTE: narrowPeak will be used by default since ATF3 is a transcription factor
module load python
pip install macs2
mkdir peakcalling
macs2 callpeak -B -t alignment/treatment1.sorted.bam -c alignment/control1.sorted.bam -f BAM -n Treatment_Control1 --nomodel --outdir peakcalling/
macs2 callpeak -B -t alignment/treatment2.sorted.bam -c alignment/control2.sorted.bam -f BAM -n Treatment_Control2 --nomodel --outdir peakcalling/

### Consensus Peaks
module load bedtools
bedtools intersect -a peakcalling/Treatment_Control1_peaks.narrowPeak -b peakcalling/Treatment_Control2_peaks.narrowPeak > peakcalling/consensus_peaks.bed

### Annotate Peaks
module load homer
awk '{$1 = "chr" $1}1' consensus_peaks.bed > consensus_peaks2.bed	# Ensures chromosome ID is in correct format
sed 's/ \+/\t/g' consensus_peaks2.bed > consensus_peaks_clean.bed	# Tab delimits file
annotatePeaks.pl consensus_peaks_clean.bed hg38  > annotations.txt

### Motif Discovery
# Extract sequences for peak coordinates
mkdir motif_discovery
bedtools getfasta -fi Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed consensus_peaks_clean.bed -s -fo motif_discovery/peak_sequences.fa

# Meme
module load perl/5.24.1
module load python/3.9.8
module load meme

meme motif_discovery/peak_sequences.fa -dna -maxw 12 -nmotifs 10 -minsites 10 -o motif_discovery/ATF3_motifs

