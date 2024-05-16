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
        usr_path = input(f'{VCSeek_config.PROMPTS[1]}')

    # Get all FASTQ files in the provided directory
    fq_list = [f for f in os.listdir(usr_path) if any(f.endswith(termination) for termination in VCSeek_config.FQ_TERMINATIONS)]

    return usr_path, fq_list, breseq


def pair_by(data_path, search_4_r1='_R1', search_4_r2='_R2'):
    """
    Pairs forward and reverse reads correctly per sample.

    Args:
        data_path (tuple): A tuple where the first element is the path and the second element is a list of filenames.
        search_4_r1 (str, optional): The substring to search for forward reads. Defaults to '_R1'.
        search_4_r2 (str, optional): The substring to replace with for reverse reads. Defaults to '_R2'.

    Returns:
        dict: A dictionary with forward reads as keys and their corresponding reverse reads as values.
    """
    
    # Initialize an empty dictionary to store paired filenames
    fq_pairs = {}
    
    # Loop through each filename in the list of filenames (data_path[1])
    for fq in data_path[1]:
        # Check if the current filename contains the forward read identifier (search_4_r1)
        if search_4_r1 in fq:
            # If the forward read identifier appears more than once in the filename, raise an error
            if fq.count(search_4_r1) > 1:
                raise NameError(f'Multiple occurrences of {search_4_r1} detected in {fq}')
            # Replace the forward read identifier with the reverse read identifier to find the pair name
            pair_name = fq.replace(search_4_r1, search_4_r2)
            # Check if the pair name exists in the list of filenames
            if pair_name in data_path[1]:
                # If it exists, add the forward and reverse read pair to the dictionary
                fq_pairs[fq] = pair_name
    
    # Return the dictionary of paired filenames
    return fq_pairs


def quali_check(data_path=None, input_dir=None, output_dir=None):
    """
    Runs FastQC on the provided data path or input directory and logs the output using MultiQC.

    This function performs quality checks on FASTQ files. It runs FastQC on each file and then
    consolidates the results using MultiQC.

    Args:
        data_path (tuple, optional): A tuple containing:
            - data_path[0] (str): The path to the directory containing the data files to be quality checked.
            - data_path[1] (list): A list of FASTQ files in the directory.
        input_dir (str, optional): The path to the input directory containing FASTQ files. Defaults to None.
        output_dir (str, optional): The path to the output directory where FastQC and MultiQC results will be saved. Defaults to None.

    Raises:
        ValueError: If neither data_path nor both input_dir and output_dir are provided.

    Example:
        Assuming `data_path` is a tuple with the first element as the directory path and the second element as a list of FASTQ files:
            data_path = ("/path/to/data", ["sample1.fastq", "sample2.fastq"])
            quali_check(data_path=data_path)

        Alternatively, you can specify `input_dir` and `output_dir`:
            quali_check(input_dir="/path/to/data", output_dir="/path/to/output")
    """
    if data_path:
        # If data_path is provided, use it to get the directory and list of FASTQ files
        data_dir = data_path[0]
        fastqc_output_dir = data_dir if output_dir is None else output_dir

        # Iterate through the list of FASTQ files
        for file in data_path[1]:
            if data_path[2]:
                breseq_fqc_out_dir = os.path.expanduser('/home/rmeneses/VSCode/Optimized_Breseq_Variant_Calling/00_quality_check_reports')
                # Construct the FastQC command
                fastqc_command = f'fastqc --noextract {data_dir}/{file} -o {breseq_fqc_out_dir}'
             
                # Run the FastQC command
                subprocess.run(fastqc_command, check=True, shell=True, capture_output=True, text=True)
              
                # Construct the MultiQC command to consolidate FastQC reports
                multiqc_command = f'multiqc {breseq_fqc_out_dir}/* -o {breseq_fqc_out_dir}'
        
            else:
                 # Construct the FastQC command
                fastqc_command = f'fastqc --noextract {data_dir}/{file} -o {fastqc_output_dir}'
            
                # Run the FastQC command
                subprocess.run(fastqc_command, check=True, shell=True, capture_output=True, text=True)

                # Construct the MultiQC command to consolidate FastQC reports
                multiqc_command = f'multiqc {fastqc_output_dir}/* -o {fastqc_output_dir}'
        
        # Run the MultiQC command
        subprocess.run(multiqc_command, check=True, shell=True, capture_output=True, text=True)

    elif input_dir and output_dir:
        # If input_dir and output_dir are provided, expand user paths
        new_input_dir = os.path.expanduser(input_dir)
        fastqc_output_dir = os.path.expanduser(output_dir)

        # Get all FASTQ files in the provided directory
        fq_list = [f for f in os.listdir(new_input_dir) if any(f.endswith(termination) for termination in VCSeek_config.FQ_TERMINATIONS)]
        print(f'FASTQ files found: {fq_list}')

        for fq_file in fq_list:
            # Construct the FastQC command
            fastqc_command = f'fastqc --noextract {new_input_dir}/{fq_file} -o {fastqc_output_dir}'
            print(f'{fastqc_command}\nProcessing: {fq_file}')
            
            # Run the FastQC command
            subprocess.run(fastqc_command, check=True, shell=True, capture_output=True, text=True)
            print(f'{fq_file} ran successfully!')

        # Construct the MultiQC command to consolidate FastQC reports
        multiqc_command = f'multiqc {fastqc_output_dir}/* -o {fastqc_output_dir}'
        print(f'Now running: {multiqc_command}')
        
        # Run the MultiQC command
        subprocess.run(multiqc_command, check=True, shell=True, capture_output=True, text=True)


def fq_pre_processing(data_path, fq_pair_dict=dict, hgt=False, rem_int=False):
    for fq_r1, fq_r2 in fq_pair_dict.items():
        sample_tag = fq_r1.split('R1')[0]
        pre_processing_output_dir = os.path.join(VCSeek_config.BASE_DIRECTORY, VCSeek_config.SUB_DIRECTORY[1])
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
        
        if hgt:
            if '_I_' in fq_r1:
                ref_genomes = VCSeek_config.ANCESTRAL_GENOMES[0]
            elif '_R_' in fq_r1:
                ref_genomes = VCSeek_config.ANCESTRAL_GENOMES[1]
        else:
            ref_genomes = f'{VCSeek_config.ANCESTRAL_GENOMES[0]},{VCSeek_config.ANCESTRAL_GENOMES[1]}'

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
        
        if rem_int:
            for int_file in intermediate_files:
                 if os.path.exists(int_file):
                    os.remove(int_file)
        print('Preprocessing completed with success!!')

        listed_processed_fq = os.listdir(pre_processing_output_dir)
    return pre_processing_output_dir, listed_processed_fq


#def alignment():

#def optmized_breseq(data_path):
    #for 