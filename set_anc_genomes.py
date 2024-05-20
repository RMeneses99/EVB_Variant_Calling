import os
from VCSeek_config import PROMPTS, ALLOWED_REF_GENOME_EXTENSIONS


def pick_ref_genome(add_ref_genome=False):
    

 if os.path.isfile(os.path.expanduser('./ref_genomes_db.py')):
        if add_ref_genome:
            ref_genome_path = input(f'{PROMPTS[2]}\n')
            while not os.path.isfile(ref_genome_path) and not any(ref_genome_path.endswith(extension) for extension in ALLOWED_REF_GENOME_EXTENSIONS):
                ref_genome_path = input(f'\n\033[31m{PROMPTS[5]}\033[0m\n')

            print(f'\n\033[0m{PROMPTS[3]}\033[0m')
            ref_genome_tag = input(f'\n\033[31m{PROMPTS[4]}\033[0m\n\n').strip()
            while not ref_genome_tag:
                ref_genome_tag = input(f'\n\033[31m{PROMPTS[6]}\033[0m\n\n').strip()

        reference_genomes[ref_genome_path] = ref_genome_tag
        
        with open('./ref_genomes_db.py', "a") as f:
            f.write(f'{reference_genomes}')
            f.close()


reference_genomes  = {}
ref_py_exists = os.path.isfile(os.path.expanduser('./ref_genomes_db.py'))
if not ref_py_exists:
    ref_genome_path = input(f'{PROMPTS[2]}\n')
    while not os.path.isfile(ref_genome_path) and not any(ref_genome_path.endswith(extension) for extension in ALLOWED_REF_GENOME_EXTENSIONS):
        ref_genome_path = input(f'\n\033[31m{PROMPTS[5]}\033[0m\n')

    print(f'\n\033[0m{PROMPTS[3]}\033[0m')
    ref_genome_tag = input(f'\n\033[31m{PROMPTS[4]}\033[0m\n\n').strip()
    while not ref_genome_tag:
        ref_genome_tag = input(f'\n\033[31m{PROMPTS[6]}\033[0m\n\n').strip()

    reference_genomes[ref_genome_path] = ref_genome_tag
        
            
    with open('./ref_genomes_db.py', "x") as f:
        f.write(f'reference_genomes = {reference_genomes }')
        f.close()

#! do you want to add another ref_genome ? Yes or No



