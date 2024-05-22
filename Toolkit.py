#!/usr/bin/python3
import os, subprocess, ast, datetime
import VCSeek_config
from ancestral_genomes_db import ancestral_genomes


#! add a feature that removes the ref folder created by bbsplit
#! ERRRORS:  hgt = input() - takes folders with no .fq/fq.gz, fastq/fastq.gz on it. You can also force your way through it without using y or n.

def ask_usr(breseq=False):
    # Initialize a list to store co-evolution information
    co_evo = []
    
    # Create the absolute path from where the script is stored
    script_path = os.path.dirname(os.path.realpath(__file__))

    # Determine the index of the base directory based on the breseq flag
    base_dir_index = 1 if breseq else 0
    
    # Determine the appropriate subdirectories based on the breseq flag
    sub_dir = VCSeek_config.BRESEQ_SUB_DIR if breseq else VCSeek_config.SUB_DIRECTORY
    
    # Create subdirectories as specified in the configuration
    for dir in sub_dir:
        os.makedirs(os.path.join(script_path,VCSeek_config.BASE_DIRECTORIES[base_dir_index], dir), exist_ok=True)

    # Prompt the user for the path containing the data to analyze
    usr_path = input(f'{VCSeek_config.PROMPTS[0]}\n')
    
    # Loop until a valid path is provided
    while not os.path.exists(usr_path):
        usr_path = input(f'\033[0m{VCSeek_config.PROMPTS[1]}\n\033[0m')

    # Ask the user if the files are from co-evolved strains
    ask_hgt = input(f'Are these files from co-evolved strains? (Y or N)\n').strip().lower()
    # If the user does not indicate a valid answer
    while ask_hgt not in ['n', 'y']:
        print('Invalid Input!!!')
        ask_hgt = input(f'Are these files from co-evolved strains? (Y or N)\n').strip().lower()

    if ask_hgt == 'y':
        tag = input(f'Please provide the ancestral reference genome tag for the isolate\n')
        while tag not in ancestral_genomes.values():
            tag = input(f'Please provide a valid ancestral reference genome tag for the isolate.\n')
        co_evo.append(tag)
        add_another_tag = input(f'Do you want to add another ancestral or plasmid in this co-evolution? (Y or N)\n').strip().lower()
        while add_another_tag not in ['n', 'y']:
            add_another_tag = input(f'Do you want to add another ancestral or plasmid in this co-evolution? (Y or N)\n').strip().lower()
        while add_another_tag == "y":
            tag = input(f'Please provide the ancestral reference genome tag for the isolate\n')
            while tag not in ancestral_genomes.values():
                tag = input(f'Please provide a valid ancestral reference genome tag for the isolate.\n')
            co_evo.append(tag)
            add_another_tag = input(f'Do you want to add another ancestral or plasmid in this co-evolution? (Y or N)\n').strip().lower()

    # Get all FASTQ files in the provided directory
    fq_list = [f for f in os.listdir(usr_path) if any(f.endswith(termination) for termination in VCSeek_config.FQ_TERMINATIONS)]
    # If for some reason duplicated tags are given, set(), method will remove them
    co_evo = set(co_evo)
    print(co_evo)
    # Return user-provided data
    return usr_path, fq_list, breseq, co_evo, script_path


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


def quali_check(data_path):
    """
    Runs FastQC on the provided data path or input directory and logs the output using MultiQC.

    This function performs quality checks on FASTQ files. It runs FastQC on each file and then
    consolidates the results using MultiQC.

    Args:
        data_path (tuple): A tuple containing:
            - data_path[0] (str): The path to the directory containing the data files to be quality checked.
            - data_path[1] (list): A list of FASTQ files in the directory.

    Raises:
        ValueError: If data_path is not provided.

    Example:
        Assuming `data_path` is a tuple with the first element as the directory path and the second element as a list of FASTQ files:
            data_path = ("/path/to/data", ["sample1.fastq", "sample2.fastq"])
            quali_check(data_path=data_path)
    """
    # If data_path is provided, use it to get the directory and list of FASTQ files
    data_dir = data_path[0]
    fastqc_output_dir = data_dir

    # Iterate through the list of FASTQ files
    for file in data_path[1]:
        # Construct the FastQC command
        fastqc_command = f'fastqc --noextract {data_dir}/{file} -o {fastqc_output_dir}'

        # Run the FastQC command
        subprocess.run(fastqc_command, check=True, shell=True, capture_output=True, text=True)

    # Construct the MultiQC command to consolidate FastQC reports
    multiqc_command = f'multiqc {fastqc_output_dir}/* -o {fastqc_output_dir}'

    # Run the MultiQC command
    subprocess.run(multiqc_command, check=True, shell=True, capture_output=True, text=True)


def pre_processing(data_path, hgt=False):
    """
    Preprocesses paired-end FASTQ files.

    Args:
        data_path (tuple): A tuple where the first element is the path and the second element is a list of filenames.
        hgt (bool, optional): Whether to perform HGT-specific processing. Defaults to False.
        remove_intermediate (bool, optional): Whether to remove intermediate files after preprocessing. Defaults to False.

    Returns:
        tuple: A tuple containing the output directory and a list of processed files.
    """

    # Use pair_by function to select pairs
    fq_pairs = pair_by(data_path)

    # Initialize variables
    processed_files = []
    if data_path[2]:
        output_directory = os.path.join(data_path[4], VCSeek_config.BASE_DIRECTORIES[1], VCSeek_config.SUB_DIRECTORY[1])
    elif not data_path[2]:
        output_directory = os.path.join(data_path[4], VCSeek_config.BASE_DIRECTORIES[0], VCSeek_config.SUB_DIRECTORY[1])

    # Iterate over each pair of forward and reverse reads
    for fq_forward, fq_reverse in fq_pairs.items():
        # Extract sample tag from forward read filename
        sample_tag = fq_forward.split('R1')[0]

        # Define fastp command
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
            '-i', f'{data_path[0]}/{fq_forward}',
            '-I', f'{data_path[0]}/{fq_reverse}',
            '-o', f'{output_directory}/processed_{fq_forward}',
            '-O', f'{output_directory}/processed_{fq_reverse}',
            '--html', '/dev/null',
            '--json', '/dev/null'
        ]
        print(f'Running: {" ".join(fastp_command)}')
        # Run fastp command
        subprocess.run(fastp_command, check=True, capture_output=True)

        # Determine reference genomes based on HGT flag
        ref_genomes = []  # Initialize ref_genomes outside of the conditional blocks
        if hgt:
        # Check if any tag from ancestral_genomes is present in the sample name
            for tag in ancestral_genomes.values():
                if tag in fq_forward:
                    ref_genomes.extend([k for k, v in ancestral_genomes.items() if v == tag])
        else:
            for tag in data_path[3]:
                ref_genomes.extend([k for k, v in ancestral_genomes.items() if v == tag])

        ref_genomes_str = ','.join(ref_genomes)

        # Define bbsplit command
        bbsplit_command = [
                          f'bbsplit.sh',
                          f'in1={output_directory}/processed_{fq_forward}',
                          f'in2={output_directory}/processed_{fq_reverse}',
                          f'ambig2=best',
                          f'ref={ref_genomes_str}',
                          f'basename={output_directory}/{sample_tag}%.fq.gz',
                          ]
        
        print(f'Running: {" ".join(bbsplit_command)}')
        # Run bbsplit command
        subprocess.run(bbsplit_command, capture_output=True)
        for i in range(len(ref_genomes)):
            ref_genome_name = ref_genomes[i].split('/')[-1]
            for ext in VCSeek_config.ALLOWED_REF_GENOME_EXTENSIONS:
                if ext in ref_genome_name:
                    ref_genome_name = ref_genome_name.split(ext)[0]
                    # Define bbsplit reformat command
                    bbsplit_reformat = ['reformat.sh',
                                       f'in={output_directory}/{sample_tag}{ref_genome_name}.fq.gz',
                                       f'out1={output_directory}/{fq_forward}',
                                       f'out2={output_directory}/{fq_reverse}'
                                       ]
#! Right now the two loops are saving in the same files: out1={output_directory}/{fq_forward} and out2={output_directory}/{fq_reverse} I need to fix this.


                    print(f'Running: {" ".join(bbsplit_reformat)}')
                    # Run bbsplit reformat command
                    subprocess.run(bbsplit_reformat, capture_output=True)

        # Print completion message
        print('Preprocessing completed successfully!')

        #! This is not sufficient has if not removed we also include the fastp output files...
        # Add processed files to list
        processed_files.extend(os.listdir(output_directory))

    # Return output directory and list of processed files
    return output_directory, processed_files

#def alignment():

#def optmized_breseq(data_path):
    #for 

def main():
    ask_usr()
    quali_check()
    pair_by()
    pre_processing()

if __name__ == '__main___':
    main()