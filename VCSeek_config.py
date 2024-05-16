import os 

ANCESTRAL_GENOMES = [os.path.expanduser("~/VSCode/Ancestral_genomes/NC_000913.gb"), os.path.expanduser("~/VSCode/Ancestral_genomes/CP054662.1_IS_annotated.gbk")]
BASE_DIRECTORIES =["./Variant_Calling", "./Optimized_Breseq_Variant_Calling"] 
SUB_DIRECTORY = ["00_quality_check_reports", "01_pre_processing", "02_alignment", "03_bam_filtering", "04_variant_calling", "05_results", "06_logs" ] 
PROMPTS = ["Please provide the path to the data you want to analyze:\n", "The provided path does not exist. Please try again.\n"]
FQ_TERMINATIONS = ['.fastq', 'fastq.gz', '.fq', 'fq.gz']
BRESEQ_SUB_DIR = ["00_quality_check_reports", "01_pre_processing", "03_breseq", "04_results", "05_logs"]