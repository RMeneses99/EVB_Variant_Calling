#!/usr/bin/python3

import os
import subprocess
import datetime
import VCSeek_config


def ask_usr(breseq=False):
    """
    Sets up necessary directories and prompts the user for a path containing FASTQ files.

    This function creates required directories based on the configuration and then prompts the user
    to input a path containing the data to analyze. It ensures the provided path is valid and returns
    the path along with a list of FASTQ files found in that directory.

    Args:
        breseq (bool, optional): If True, creates directories specific to Breseq analysis. Defaults to False.

    Returns:
        tuple: A tuple containing:
            - usr_path (str): The user-provided path containing the data files.
            - fq_list (list): A list of FASTQ files found in the provided directory.
            - breseq (bool): The boolean value indicating if Breseq-specific directories were created.

    Raises:
        ValueError: If the user provides an invalid path.
    """
    if breseq:
        # If breseq flag is True, create Breseq-specific subdirectories
        for directory in VCSeek_config.BRESEQ_SUB_DIR:
            os.makedirs(os.path.join(VCSeek_config.BASE_DIRECTORIES[1], directory), exist_ok=True)
    else:
        # Create general subdirectories as specified in the configuration
        for directory in VCSeek_config.SUB_DIRECTORY:
            # Create the directory if it does not exist already
            os.makedirs(os.path.join(VCSeek_config.BASE_DIRECTORIES[0], directory), exist_ok=True)

    # Prompt the user for the path containing the data to analyze
    usr_path = input(f'{VCSeek_config.PROMPTS[0]}')

    # Loop until a valid path is provided
    while not os.path.exists(usr_path):
        usr_path = input(f'\033[0m{VCSeek_config.PROMPTS[1]}\033[0m')

    # Get all FASTQ files in the provided directory
    fq_list = [f for f in os.listdir(usr_path) if any(f.endswith(termination) for termination in VCSeek_config.FQ_TERMINATIONS)]

    return usr_path, fq_list, breseq

    
def quali_check(data_path=None, input_dir=None, output_dir=None, search_4_r1='_R1', search_4_r2='_R2'):
    # Initialize an empty dictionary to store paired FASTQ files
    fq_pairs = {}
    
    # If data_path is provided
    if data_path:
        # Loop through the list of FASTQ files in data_path
        for fq in data_path[1]:
            # Check if the file contains the search pattern for the first read
            if search_4_r1 in fq:
                # Raise an error if multiple occurrences of the pattern are found
                if fq.count(search_4_r1) > 1:
                    raise NameError(f'Multiple occurrences of {search_4_r1} detected in {fq}')
            
            # Replace the pattern for the first read with the pattern for the second read to find the pair
            pair_name = fq.replace(search_4_r1, search_4_r2)
            
            # Check if the pair exists in the list of files and add to the dictionary
            if pair_name in data_path[1]:
                fq_pairs[fq] = pair_name

        # Set the data directory and output directory for FastQC
        data_dir = data_path[0]
        fastqc_output_dir = data_dir if output_dir is None else output_dir
        
        # Loop through each file in data_path
        for file in data_path[1]:
            if data_path[2]:
                # Use specific output directory for FastQC reports if indicated
                breseq_fqc_out_dir = os.path.expanduser('/home/rmeneses/VSCode/Optimized_Breseq_Variant_Calling/00_quality_check_reports')
                fastqc_command = f'fastqc --noextract {data_dir}/{file} -o {breseq_fqc_out_dir}'
                subprocess.run(fastqc_command, check=True, shell=True, capture_output=True, text=True)
                multiqc_command = f'multiqc {breseq_fqc_out_dir}/* -o {breseq_fqc_out_dir}'
            else:
                # Use the provided or default output directory for FastQC reports
                fastqc_command = f'fastqc --noextract {data_dir}/{file} -o {fastqc_output_dir}'
                subprocess.run(fastqc_command, check=True, shell=True, capture_output=True, text=True)
                multiqc_command = f'multiqc {fastqc_output_dir}/* -o {fastqc_output_dir}'
        
        # Run MultiQC to aggregate FastQC reports
        subprocess.run(multiqc_command, check=True, shell=True, capture_output=True, text=True)
    
    # If input_dir and output_dir are provided
    elif input_dir and output_dir:

        # Expand user directory for input and output directories
        new_input_dir = os.path.expanduser(input_dir)
        fastqc_output_dir = os.path.expanduser(output_dir)
        
        # Ensure input directory exists
        if not os.path.isdir(new_input_dir):
            raise FileNotFoundError(f"Input directory {new_input_dir} does not exist")
        
        # Ensure output directory exists, create if not
        if not os.path.isdir(fastqc_output_dir):
            os.makedirs(fastqc_output_dir, exist_ok=True)
        
        # List FASTQ files in the input directory based on specific terminations
        fq_list = [f for f in os.listdir(new_input_dir) if any(f.endswith(termination) for termination in VCSeek_config.FQ_TERMINATIONS)]
        print(f'FASTQ files found: {fq_list}')
        
        # Run FastQC for each FASTQ file
        for fq_file in fq_list:
            fastqc_command = f'fastqc --noextract {new_input_dir}/{fq_file} -o {fastqc_output_dir}'
            print(f'{fastqc_command}\nProcessing: {fq_file}')
            subprocess.run(fastqc_command, check=True, shell=True, capture_output=True, text=True)
            print(f'{fq_file} ran successfully!')
        
        # Run MultiQC to aggregate FastQC reports
        multiqc_command = f'multiqc {fastqc_output_dir}/* -o {fastqc_output_dir}'
        print(f'Now running: {multiqc_command}')
        subprocess.run(multiqc_command, check=True, shell=True, capture_output=True, text=True)
    else:
        print(f'')


def fq_pre_processing(data_path=None, hgt=False, rm_inter_files=False, input_dir=None, output_dir=None, search_4_r1='_R1', search_4_r2='_R2'):
    #! hgt=False is counter intuitive
    fq_pairs = {}

    if data_path:
       # Loop through the list of FASTQ files in data_path
        for fq in data_path[1]:
            # Check if the file contains the search pattern for the first read
            if search_4_r1 in fq:
                # Raise an error if multiple occurrences of the pattern are found
                if fq.count(search_4_r1) > 1:
                    raise NameError(f'Multiple occurrences of {search_4_r1} detected in {fq}')
                # Replace the pattern for the first read with the pattern for the second read to find the pair
                pair_name = fq.replace(search_4_r1, search_4_r2)
                # Check if the pair exists in the list of files and add to the dictionary
                if pair_name in data_path[1]:
                    fq_pairs[fq] = pair_name
            for fq_r1, fq_r2 in fq_pairs.list():
                sample_tag = fq_r1.split('R1')[0]
                pre_processing_output_dir = os.path.join(VCSeek_config.BASE_DIRECTORY, VCSeek_config.SUB_DIRECTORY[1])
                if hgt:
                    if '_I_' in fq_r1:
                        ref_genomes = VCSeek_config.ANCESTRAL_GENOMES[0]
                    elif '_R_' in fq_r1:
                        ref_genomes = VCSeek_config.ANCESTRAL_GENOMES[1]
                    else:
                        print(f'Your sample does not have _I_ or _R_ identifiers, please add them to your sample name.')
                else:
                    ref_genomes = f'{VCSeek_config.ANCESTRAL_GENOMES[0]},{VCSeek_config.ANCESTRAL_GENOMES[1]}'
                fastp_command = [
                    'fastp',
                    '-q', '20',
                    '-u', '50',
                    '--dedup', '1',
                    '--detect_adapter_for_pe',
                    '-p', '3',
                    '-5',
                    '-M', '20',
                    '-W', '4',
                    '-c',
                    '-i', f'{data_path[0]}/{fq_r1}',
                    '-I', f'{data_path[0]}/{fq_r2}',
                    '-o', f'{pre_processing_output_dir}/processed_{fq_r1}',
                    '-O', f'{pre_processing_output_dir}/processed_{fq_r2}',
                    '--html', '/dev/null',
                    '--json', '/dev/null'
                    ]

                subprocess.run(fastp_command, check=True)

                bbsplit_command = [
                'bbsplit.sh',
                f'in1={pre_processing_output_dir}/processed_{fq_r1}',
                f'in2={pre_processing_output_dir}/processed_{fq_r2}',
                'ambig2=best',
                f'ref={ref_genomes}',
                f'basename={pre_processing_output_dir}/{sample_tag}%.fq.gz'
                ]   
    
                subprocess.run(bbsplit_command, check=True, capture_output=True, text=True)
                bbsplit_reformat = ['reformat.sh',
                f'in={pre_processing_output_dir}/{sample_tag}NC_000913.fq.gz',
                f'out1={pre_processing_output_dir}/{fq_r1}',
                f'out2={pre_processing_output_dir}/{fq_r2}'
                ]
                subprocess.run(bbsplit_reformat, check=True, capture_output=True)

                intermediate_files = [f'{pre_processing_output_dir}/processed_{fq_r1}', 
                                      f'{pre_processing_output_dir}/processed_{fq_r2}', 
                                      f'{pre_processing_output_dir}/{sample_tag}NC_000913.fq.gz' 
                                      ]
        
                if rm_inter_files:
                    for int_file in intermediate_files:
                        if os.path.exists(int_file):
                            os.remove(int_file)
        print('Preprocessing completed with success!!')

    elif input_dir and output_dir:
        for fq in data_path[1]:
            # Check if the file contains the search pattern for the first read
            if search_4_r1 in fq:
                # Raise an error if multiple occurrences of the pattern are found
                if fq.count(search_4_r1) > 1:
                    raise NameError(f'Multiple occurrences of {search_4_r1} detected in {fq}')
                # Replace the pattern for the first read with the pattern for the second read to find the pair
                pair_name = fq.replace(search_4_r1, search_4_r2)
                # Check if the pair exists in the list of files and add to the dictionary
                if pair_name in data_path[1]:
                    fq_pairs[fq] = pair_name
            for fq_r1, fq_r2 in fq_pairs.list():
                sample_tag = fq_r1.split('R1')[0]
                pre_processing_output_dir = os.path.join(VCSeek_config.BASE_DIRECTORY, VCSeek_config.SUB_DIRECTORY[1])
                if hgt:
                    if '_I_' in fq_r1:
                        ref_genomes = VCSeek_config.ANCESTRAL_GENOMES[0]
                    elif '_R_' in fq_r1:
                        ref_genomes = VCSeek_config.ANCESTRAL_GENOMES[1]
                    else:
                        ref_genomes = f'{VCSeek_config.ANCESTRAL_GENOMES[0]},{VCSeek_config.ANCESTRAL_GENOMES[1]}'
                fastp_command = [
                    'fastp',
                    '-q', '20',
                    '-u', '50',
                    '--dedup', '1',
                    '--detect_adapter_for_pe',
                    '-p', '3',
                    '-5',
                    '-M', '20',
                    '-W', '4',
                    '-c',
                    '-i', f'{data_path[0]}/{fq_r1}',
                    '-I', f'{data_path[0]}/{fq_r2}',
                    '-o', f'{pre_processing_output_dir}/processed_{fq_r1}',
                    '-O', f'{pre_processing_output_dir}/processed_{fq_r2}',
                    '--html', '/dev/null',
                    '--json', '/dev/null'
                    ]

                subprocess.run(fastp_command, check=True)

                bbsplit_command = [
                'bbsplit.sh',
                f'in1={pre_processing_output_dir}/processed_{fq_r1}',
                f'in2={pre_processing_output_dir}/processed_{fq_r2}',
                'ambig2=best',
                f'ref={ref_genomes}',
                f'basename={pre_processing_output_dir}/{sample_tag}%.fq.gz'
                ]   
    
                subprocess.run(bbsplit_command, check=True, capture_output=True, text=True)
                bbsplit_reformat = ['reformat.sh',
                f'in={pre_processing_output_dir}/{sample_tag}NC_000913.fq.gz',
                f'out1={pre_processing_output_dir}/{fq_r1}',
                f'out2={pre_processing_output_dir}/{fq_r2}'
                ]
                subprocess.run(bbsplit_reformat, check=True, capture_output=True)

                intermediate_files = [f'{pre_processing_output_dir}/processed_{fq_r1}', 
                                      f'{pre_processing_output_dir}/processed_{fq_r2}', 
                                      f'{pre_processing_output_dir}/{sample_tag}NC_000913.fq.gz' 
                                      ]
        
                if rm_inter_files:
                    for int_file in intermediate_files:
                        if os.path.exists(int_file):
                            os.remove(int_file)
        print('Preprocessing completed with success!!')
    
    listed_processed_fq = os.listdir(pre_processing_output_dir)
    return pre_processing_output_dir, listed_processed_fq


#def alignment():

#def optmized_breseq(data_path):
    #for 