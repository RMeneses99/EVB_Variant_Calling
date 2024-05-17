import os 

#ANCESTRAL_GENOMES = [os.path.expanduser("~/VSCode/Ancestral_genomes/NC_000913.gb"), os.path.expanduser("~/VSCode/Ancestral_genomes/CP054662.1_IS_annotated.gbk")]
ANCESTRAL_GENOMES = {
                    os.path.expanduser("~/VSCode/Ancestral_genomes/NC_000913.gb"):'_I_',
                    os.path.expanduser("~/VSCode/Ancestral_genomes/CP054662.1_IS_annotated.gbk") : '_R_'  
                    }

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
          "Note: that the reference genome identifier tag should be present somewhere in the sample name.",
          "The provided path is not for a file or does not have an accepted extension. Please provide a valid file",
          "Your reference genome identification tag cannot be null or an empty space"
          ]

FQ_TERMINATIONS = [
                  '.fastq',
                  'fastq.gz',
                  '.fq',
                  'fq.gz'
                  ]

ALLOWED_REF_GENOME_EXTENSIONS = [".fa",".fasta",".gff3", ".gbk", ".gk"]

BRESEQ_SUB_DIR = [
                 "00_quality_check_reports",
                 "01_pre_processing",
                 "03_breseq",
                 "04_results",
                 "05_logs"
                 ]