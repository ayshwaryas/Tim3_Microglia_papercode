# tim3_microglia
## Bulk-RNA seq
### 60 samples
### GET site data(68 samples)
Data files:
 - Merged raw fastq files - /broad/kuchroolab/kimi_microglia/GET_data/get.broadinstitute.org/pkgs/SN0238404/merged_fastqs/
 - RSEM result files- /broad/kuchroolab/kimi_microglia/GET_data/rsem_data/
 - Sample sheets and metadata- /broad/kuchroolab/kimi_microglia/GET_data/samplesheet+metadata/
 - Count and TPM metrices .csv files- /broad/kuchroolab/kimi_microglia/GET_data/


## Single nucleus (snRNA-seq)

Data files:
- Raw 10X data- /broad/kuchroolab/kimi_microglia/nucseq/215444255/
- Sample sheet for pre-processing - /broad/kuchroolab/kimi_microglia/nucseq/pre-processing/nucseq_samplesheet_kimi_microglia.csv
- Cellranger 10X output- /broad/kuchroolab/kimi_microglia/nucseq/raw_data/

Pre-processing- cellranger_workflow against pre-mRNA( mm10_premrna) reference genome.
Downstream Analysis- QC, Harmony integration, clustering, sub-clustering, DEGs per cluster, and annotations


