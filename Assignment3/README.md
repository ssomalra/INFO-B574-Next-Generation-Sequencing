# Assignment 3
1. Run the ChIP-seq analysis on the demo data from H3K4me3 and Control:
   
   [`H3K4m31_R1`](http://tracks.ccbb.iupui.edu/tracks/sliu19/NGS2024/H3K4me3_S1_R1_L001.fastq.gz)<br>
   [`H3K4m31_R2`](http://tracks.ccbb.iupui.edu/tracks/sliu19/NGS2024/H3K4me3_S1_R2_L001.fastq.gz)<br>
   [`Control_R1`](http://tracks.ccbb.iupui.edu/tracks/sliu19/NGS2024/Control_R1_L001.fastq.gz)<br>
   [`Control_R2`](http://tracks.ccbb.iupui.edu/tracks/sliu19/NGS2024/Control_R2_L001.fastq.gz)
2. Generate the fragment size distribution figure and fingerprint figure for both samples
3. Call peaks using MACS3 with two options, one narrowPeak and one broadPeak. Generate the xls files and tell the difference between the two.

## Assignment Submission
**Fragment size distribution figures**
<div style="display: flex; justify-content: space-between;">
    <img src="https://github.com/user-attachments/assets/0ad15026-2e5d-4476-995d-93fdfd6a0dc3" alt="Control_fragment" width="45%">
    <img src="https://github.com/user-attachments/assets/17374c95-94a0-4aa8-9fa0-5f304733f74c" alt="H3K4me3_S1_fragment" width="45%">
</div>

**Fingerprint figures**
<div style="display: flex; justify-content: space-between;">
    <img src="https://github.com/user-attachments/assets/1927cde3-597d-4ba6-8cf9-b5ac9f9493a7" alt="Control_fingerprint" width="45%">
    <img src="https://github.com/user-attachments/assets/db69b83e-ea33-40cc-9fc8-0801d5172d7b" alt="H3K4me3_S1_fingerprint" width="45%">
</div>

<br>**narrowPeak vs broadPeak**

The primary distinction lies in the extent of the peak being captured; specifically, the narrowPeak file encompasses a more restricted region bound compared to the broadPeak file.

Additionally, the narrowPeak file contains the "abs_summit" feature, which represents the precise genomic coordinate where the signal intensity reaches its maximum within a peak. This metric is particularly useful for pinpointing the exact locations of regulatory elements like transcription factor binding sites. In contrast, the broadPeak file lacks this attribute as it represents broader regions of enrichment rather than precise peaks.

Finally, the narrowPeak file typically exhibits higher fold enrichment scores compared to broadPeak files. This discrepancy may stem from narrowPeak regions capturing specific binding events with concentrated signals within a narrow genomic region. This focused signal results in higher read counts within the peak region compared to the background, leading to higher fold enrichment values compared to the broadPeak regions, where signal distribution tends to be more diffuse across broader genomic regions.
