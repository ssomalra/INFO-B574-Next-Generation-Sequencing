##GATK short variants discovery pipeline

##get the demo dataset
mkdir -p /N/scratch/$USER/GATK4
cd /N/scratch/$USER/GATK4

wget --recursive --continue --no-host-directories --no-parent --cut-dirs=1 \
    -A "*.gz,*.bam,*.bai,*.json,*.txt,*.html,*.h5,*.cloupe,*tbi,*.csv" \
    -P G788_spring2024 https://cmg-data.sca.iu.edu/2d8cf0e6-8246-4305-a966-6d8207f6f2f2

tar xvf G788_spring2024/WGS_data.tar.gz

#================================================================================
##load modules
module load fastqc
module load bwa
module load samtools
module load java/15.0.2 gatk
module load picard
picard=/N/soft/rhel8/picard/3.0.0/build/libs/picard.jar

#================================================================================
#QC
cd /N/scratch/$USER/GATK4

#FASTQC
module load fastqc
mkdir -p report/fastqc/rawfq
fastqc FASTQ/*.fastq.gz -o report/fastqc/rawfq/

# ADAPTER TRIMMING
mkdir -p FASTQ/trimming/
module load trimmomatic
# for sample 1
trimmomatic PE -phred33 \
    FASTQ/S1_R1_001.fastq.gz FASTQ/S1_R2_001.fastq.gz \
    FASTQ/trimming/S1_R1_001_paired.fq.gz FASTQ/trimming/S1_R1_001_unpaired.fq.gz \
    FASTQ/trimming/S1_R2_001_paired.fq.gz FASTQ/trimming/S1_R2_001_unpaired.fq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:50
# for sample 2    
trimmomatic PE -phred33 \
    FASTQ/S2_R1_001.fastq.gz FASTQ/S2_R2_001.fastq.gz \
    FASTQ/trimming/S2_R1_001_paired.fq.gz FASTQ/trimming/S2_R1_001_unpaired.fq.gz \
    FASTQ/trimming/S2_R2_001_paired.fq.gz FASTQ/trimming/S2_R2_001_unpaired.fq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:50

#FASTQC
module load fastqc
mkdir -p report/fastqc/trimfq
fastqc FASTQ/trimming/*_paired.fq.gz -o report/fastqc/trimfq/

#================================================================================
#ALIGNMENT
GATK_bundle=Reference/broad.hg38.v0/
reference=$GATK_bundle/Homo_sapiens_assembly38.fasta
NT=8
SID=S1 # first run alignment for the first sample
# SID=S2 # then run alignment for the second sample
read_group="@RG\tID:${SID}\tSM:${SID}\tPL:ILLUMINA"

mkdir -p report/alignment/${SID}/
bwa mem \
    -t $NT \
    -R ${read_group} \
    ${reference} \
    FASTQ/trimming/${SID}_R1_001_paired.fq.gz \
    FASTQ/trimming/${SID}_R2_001_paired.fq.gz \
    > report/alignment/${SID}/${SID}.sam

#SORTING AND INDEXING
samtools view -bS report/alignment/${SID}/${SID}.sam > report/alignment/${SID}/${SID}.bam 
samtools sort -@ $NT -o report/aligmnment/${SID}/${SID}.sorted.bam report/alignment/${SID}/${SID}.bam 
samtools index report/alignment/${SID}/${SID}.sorted.bam

#MARK DUPLICATES
module load picard
picard=/N/soft/rhel8/picard/3.0.0/build/libs/picard.jar
java -jar $picard MarkDuplicates \
    I=report/alignment/${SID}/${SID}.sorted.bam \
    O=report/alignment/${SID}/${SID}.sorted.mkdp.bam \
    M=report/alignment/${SID}/${SID}.sorted.mkdp_metrics.txt

samtools index report/alignment/${SID}/${SID}.sorted.mkdp.bam

#================================================================================
#BQSR
module load gatk
SID=S1
# SID=S2

gatk BaseRecalibrator \
    -R ${reference} \
    -I report/alignment/${SID}/${SID}.sorted.mkdp.bam \
    --known-sites $GATK_bundle/Homo_sapiens_assembly38.dbsnp138.vcf \
    --known-sites $GATK_bundle/Homo_sapiens_assembly38.known_indels.vcf.gz \
    -O report/alignment/${SID}/${SID}.recal_data.table
 
gatk ApplyBQSR \
   -R ${reference} \
   -I report/alignment/${SID}/${SID}.sorted.mkdp.bam \
   --bqsr-recal-file report/alignment/${SID}/${SID}.recal_data.table \
   -O report/alignment/${SID}/${SID}.sorted.mkdp.bqsr.bam

#================================================================================

#variant calling for one-sample analysis
module load gatk
SID=S1
# SID=S2

mkdir -p report/gatk/${SID}

# Run HaplotypeCaller on a single sample
gatk HaplotypeCaller \
  -R ${reference} \
  -I report/alignment/${SID}/${SID}.sorted.mkdp.bqsr.bam \
  -ERC GVCF \
  -O report/gatk/${SID}/${SID}.g.vcf.gz

    #plan A: if there is only one sample, run GenotypeGVCFs on a single sample
    gatk GenotypeGVCFs \
      -R ${reference} \
      -V report/gatk/${SID}/${SID}.g.vcf.gz \
      -O report/gatk/${SID}/${SID}.vcf.gz

    # Apply hard filtering
    gatk VariantFiltration \
      -R ${reference} \
      -V report/gatk/${SID}/${SID}.vcf.gz \
      --filter 'vc.isSNP() && (QD<2.0 || FS>60.0 || MQ<40.0)' \
      --filter-name "SNP_FILTER" \
      --filter '!vc.isSNP() && (QD<2.0 || FS>200.0 || MQ<30.0)' \
      --filter-name "INDEL_FILTER" \
      -O report/gatk/${SID}/${SID}.hard.vcf.gz

#================================================================================

##variant calling for large scale study

##Run alignment and BQSR for all samples
##Run the code above before "HaplotypeCaller" with SID=S2

##Run HaplotypeCaller on each sample
for SID in S1 S2
do
    if [[ ! -f report/gatk/${SID}/${SID}.g.vcf.gz ]];then 
        gatk HaplotypeCaller \
          -R ${reference} \
          -I report/alignment/${SID}/${SID}.sorted.mkdp.bqsr.bam \
          -ERC GVCF \
          -O report/gatk/${SID}/${SID}.g.vcf.gz
    fi
done

##run by chromosome separately. this dataset has only chr11 reads
test -e report/gatk/GenomicsDB && rm -rf report/gatk/GenomicsDB #delete the folder if it exists
gatk GenomicsDBImport \
    -V report/gatk/S1/S1.g.vcf.gz \
    -V report/gatk/S2/S2.g.vcf.gz \
    -L chr11 \
    --genomicsdb-workspace-path report/gatk/GenomicsDB

gatk GenotypeGVCFs \
    -R ${reference} \
    -V gendb://report/gatk/GenomicsDB \
    -O report/gatk/joint_genotyped.vcf.gz

##Run GATK VariantRecalibrator for SNPs
input_vcf=report/gatk/joint_genotyped.vcf.gz
mkdir -p report/gatk/VQSR

gatk VariantRecalibrator \
    -R $reference \
    -V $input_vcf \
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $GATK_bundle/hapmap_3.3.hg38.vcf.gz \
    -resource:omni,known=false,training=true,truth=false,prior=12.0 $GATK_bundle/1000G_omni2.5.hg38.vcf.gz \
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 $GATK_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $GATK_bundle/Homo_sapiens_assembly38.dbsnp138.vcf \
    -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum \
    -mode SNP \
    -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
    --tranches-file report/gatk/VQSR/snp.tranches \
    --rscript-file report/gatk/VQSR/output_SNP.plots.R \
    -O report/gatk/VQSR/snp.recal

##Apply VQSR to the SNPs in the VCF
gatk ApplyVQSR \
    -R $reference \
    -V $input_vcf \
    --recal-file report/gatk/VQSR/snp.recal \
    --tranches-file report/gatk/VQSR/snp.tranches \
    --truth-sensitivity-filter-level 99.0 \
    -mode SNP \
    -O report/gatk/VQSR/snp.output.vcf

##Run GATK VariantRecalibrator for INDELs
gatk VariantRecalibrator \
    -R $reference \
    -V report/gatk/VQSR/snp.output.vcf \
    -resource:mills,known=false,training=true,truth=true,prior=12.0 $GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $GATK_bundle/Homo_sapiens_assembly38.known_indels.vcf.gz \
    -an DP -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum \
    -mode INDEL \
    -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
    --tranches-file report/gatk/VQSR/indel.tranches \
    --rscript-file report/gatk/VQSR/output_indel.plots.R \
    -O report/gatk/VQSR/indel.recal

##Apply VQSR to the INDELs in the VCF
gatk ApplyVQSR \
    -R $reference \
    -V report/gatk/VQSR/snp.output.vcf \
    --recal-file report/gatk/VQSR/indel.recal \
    --tranches-file report/gatk/VQSR/indel.tranches \
    --truth-sensitivity-filter-level 99.0 \
    -mode INDEL \
    -O report/gatk/final.output.vcf

#================================================================================

##variant annotation using ANNOVAR

##for one-sample analysis
gunzip report/gatk/S1/S1.hard.vcf.gz -c > report/gatk/S1/S1.hard.vcf
input_vcf=report/gatk/S1/S1.hard.vcf

##for multi-sample analysis
# input_vcf=report/gatk/final.output.vcf

mkdir report/ANNOVAR

table_annovar.pl \
  $input_vcf \
  hg38/ \
  -buildver hg38 \
  -out report/ANNOVAR/WGS \
  -remove \
  -protocol refGene,cytoBand,exac03,avsnp150,dbnsfp41a \
  -operation g,r,f,f,f \
  -nastring . \
  -vcfinput

