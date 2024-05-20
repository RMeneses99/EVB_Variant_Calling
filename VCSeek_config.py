import os 


BASE_DIRECTORIES = [
                   "./Variant_Calling", 
                   "./Optimized_Breseq_Variant_Calling"
                   ] 

SUB_DIRECTORY = [
                "00_quality_check_reports", 
                "01_pre_processing", 
                "02_alignment", 
                "03_bam_filtering", 
                "04_variant_calling", 
                "05_results", 
                "06_logs" 
                ] 

PROMPTS = [
          "Please provide the path to the data you want to analyze:",
          "The provided path does not exist. Please try again:",
          "Please provide the path to the genome file you want to use has reference:",
          "Please provide the identifier you want to attribute to this reference genome.",
          "Note: The reference genome identifier tag should be somwhere present in the sample name.",
          "The provided path is not for a file or does not have an accepted extension. Please provide a valid file",
          "Your reference genome identification tag cannot be null or an empty space"
          ]

FQ_TERMINATIONS = [
                  '.fastq',
                  'fastq.gz',
                  '.fq',
                  'fq.gz'
                  ]

ALLOWED_REF_GENOME_EXTENSIONS = [
                                ".fa",
                                ".fasta",
                                ".gff3",
                                ".gbk",
                                ".gb"]

BRESEQ_SUB_DIR = [
                 "00_quality_check_reports", 
                 "01_pre_processing", 
                 "02_alignment", 
                 "03_bam_filtering",
                 "05_results", 
                 "06_logs" 
                 ]

TEST = {}