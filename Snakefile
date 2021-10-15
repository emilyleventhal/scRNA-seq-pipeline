
GEO_SAMPLES = ["SRR13694201", "SRR13694213","SRR13694212", "SRR13694214", 
"SRR13694202", "SRR13694200", "SRR13694205", "SRR13694204", "SRR13694206", \
"SRR13694207", "SRR13694208", "SRR13694210", "SRR13694211", "SRR13694209",  \
"SRR13694215"] 

FILES = ["MS3686b1_JH_L1-ds.13d32cd0631a4142b34e1b79ee52dd98",  "MS3686b3_JH_L1-ds.1a6c2802dd0e4cec93f899c94eb00670", "MS3686b2_JH_L1-ds.eecbf39e34ef48d3a4121ff147d20ec3", "MS3686b4_JH_L1-ds.668ecbf5bc5140ecb738f12a13641c54"]

SAMPLES = ["sample1", "sample2", "sample3", "sample4"]


rule all: 
    input: expand("GEO_datasets/{geo_sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", geo_sample = GEO_SAMPLES), \
    expand("/scratch/ell2yvj/projects/scRNA-seq/GEO_datasets/{geo_sample}/{geo_sample}_S1_L001_R1_001_fastqc.html", geo_sample = GEO_SAMPLES)

rule get_fastq: 
    output: "GEO_datasets/{geo_sample}/{geo_sample}.fastq"
    threads: 1 
    resources: mem_mb = 30000
    shell: """
    export PATH=$PATH:/scratch/ell2yvj/projects/scATAC-seq/sratoolkit.2.11.1-centos_linux64/bin ; 
    fastq-dump {wildcards.geo_sample} -O /scratch/ell2yvj/projects/scRNA-seq/GEO_datasets/{wildcards.geo_sample}/{wildcards.geo_sample}.fastq ;
    rm /home/ell2yvj/user-repo/sra/{wildcards.geo_sample}.sra
    """

rule split_fastq: 
    input: "GEO_datasets/{geo_sample}/{geo_sample}.fastq"
    output: "GEO_datasets/{geo_sample}/{geo_sample}_1.fastq.gz", \
    "GEO_datasets/{geo_sample}/{geo_sample}_2.fastq.gz", \
    "GEO_datasets/{geo_sample}/{geo_sample}_3.fastq.gz"
    threads: 1
    resources: mem_mb=30000
    shell: """cd /scratch/ell2yvj/projects/scRNA-seq/GEO_datasets/{wildcards.geo_sample} ;
    export PATH=$PATH:/scratch/ell2yvj/projects/scATAC-seq/sratoolkit.2.11.1-centos_linux64/bin ;
    fastq-dump --split-files --gzip {wildcards.geo_sample} 
    """

rule fastqc: 
   threads: 1 
   resources: mem_mb=50000
   input: "GEO_datasets/{geo_sample}/{geo_sample}_S1_L001_I1_001.fastq.gz", \
   "GEO_datasets/{geo_sample}/{geo_sample}_S1_L001_R1_001.fastq.gz", \
   "GEO_datasets/{geo_sample}/{geo_sample}_S1_L001_R2_001.fastq.gz"
   output: "GEO_datasets/{geo_sample}/{geo_sample}_S1_L001_I1_001_fastqc.html", \
   "GEO_datasets/{geo_sample}/{geo_sample}_S1_L001_R1_001_fastqc.html", \
   "GEO_datasets/{geo_sample}/{geo_sample}_S1_L001_R2_001_fastqc.html"
   shell: """ 
   module load fastqc ; 
   fastqc /scratch/ell2yvj/projects/scRNA-seq/GEO_datasets/{wildcards.geo_sample}/*.fastq.gz -o /scratch/ell2yvj/projects/scRNA-seq/GEO_datasets/{wildcards.geo_sample}
   """

# need three files with these exact names, and you must unzip them: 
# barcodes.tsv
# genes.tsv
# matrix.mtx

## named w/ wrong R1, R2, I1 
rule change_name_1:
    input: "GEO_datasets/{geo_sample}/{geo_sample}_1.fastq.gz" # change to {geo_sample}_1.fastq.gz eventually
    output: "GEO_datasets/{geo_sample}/{geo_sample}_S1_L001_I1_001.fastq.gz"
    shell: """
    mv {input} {output}
    """
rule change_name_2:
    input: "GEO_datasets/{geo_sample}/{geo_sample}_2.fastq.gz" # change to {geo_sample}_2.fastq.gz
    output: "GEO_datasets/{geo_sample}/{geo_sample}_S1_L001_R1_001.fastq.gz"
    shell: """
    mv {input} {output}
    """
     #SampleName_S1_L001_R1_001.fastq.gz
rule change_name_3: 
    input: "GEO_datasets/{geo_sample}/{geo_sample}_3.fastq.gz"
    output: "GEO_datasets/{geo_sample}/{geo_sample}_S1_L001_R2_001.fastq.gz"
    shell: """
    mv {input} {output}
    """

rule count: 
    input: "GEO_datasets/{geo_sample}/{geo_sample}_S1_L001_I1_001.fastq.gz", \
    "GEO_datasets/{geo_sample}/{geo_sample}_S1_L001_R1_001.fastq.gz", \
    "GEO_datasets/{geo_sample}/{geo_sample}_S1_L001_R2_001.fastq.gz"
    output: "/scratch/ell2yvj/projects/scRNA-seq/{geo_sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
    threads: 1
    resources: mem_mb = 64000
    shell: """
        export PATH=/scratch/ell2yvj/projects/scRNA-seq/cellranger-6.1.1:$PATH ;
        rm -r {wildcards.geo_sample} ; 
        cellranger count --id={wildcards.geo_sample} \
        --fastqs=/scratch/ell2yvj/projects/scRNA-seq/GEO_datasets/{wildcards.geo_sample}/ \
        --transcriptome=/project/ScottLab/Mouse.ZsG.genome/ \
        --include-introns \
        --localcores=16
        """

# rule count_our_samples: 
#     output: "/scratch/ell2yvj/projects/scRNA-seq/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
#     threads: 1
#     resources: mem_mb = 64000
#     shell: """
#         export PATH=/scratch/ell2yvj/projects/scRNA-seq/cellranger-6.1.1:$PATH ;
#         rm -r {sample} ; 
#         cellranger count --id={wildcards.sample} \
#         --fastqs=../JH_scRNASeq/M089ScottMS3686DudleyYZ364810x3p-456194756/{wildcards.sample}/ \
#         --transcriptome=/project/ScottLab/Mouse.ZsG.genome/ \
#         --include-introns \
#         --localcores=16
#         """
        # export PATH=/scratch/ell2yvj/projects/scRNA-seq/cellranger-6.1.1:$PATH
        # cellranger count --id=fulldata_S1\
        #     --fastqs=/scratch/ell2yvj/projects/scRNA-seq/rawdata/JH_scRNASeq/FASTQ_Generation_2021-08-31_00_45_57Z-455128675/M088MS3686_01_L001-ds.1e6a2020d0574d35be5ed8852b18b218 \
        #     --transcriptome=/scratch/ell2yvj/projects/scRNA-seq/Mus_musculus.GRCm38.98.premrna.filtered
        # """
        
# rule getPlots: 
#     input: data_dir = "rawdata/GSM3559978/"
#     output: QC_metrics = "outputs/GSM3559978/QC_metrics.txt", \
#     QC_metrics_vln = "outputs/GSM3559978/QC_metrics_vln.pdf", \
#     feature_feature = "outputs/GSM3559978/feature-feature-relationships.pdf", \
#     top_var = "outputs/GSM3559978/top_variable_genes.txt", \
#     var_features = "outputs/GSM3559978/variable_features.pdf", \
#     scaled_data = "outputs/GSM3559978/scaled_data.txt", \
#     PCA_results = "outputs/GSM3559978/PCA_results.txt"
#     "outputs/GSM3559978/PCA_visualizations"
