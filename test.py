import os
import ast
from VCSeek_config import PROMPTS, ALLOWED_REF_GENOME_EXTENSIONS


def pick_ref_genome(add_ref_genome=False):

    ref_genomes_py = os.path.expanduser('./ref_genomes_db.py')
    ref_py_exists = os.path.isfile(ref_genomes_py)
    reference_genomes = {}

    # Load existing reference genomes if the file exists
    if ref_py_exists:
        with open(ref_genomes_py, 'r') as f:
            content = f.read()
            tree = ast.parse(content, mode='exec')
            for node in tree.body:
                if isinstance(node, ast.Assign):
                    for target in node.targets:
                        if isinstance(target, ast.Name) and target.id == "reference_genomes":
                            reference_genomes = ast.literal_eval(node.value)

    while add_ref_genome:
        ref_genome_path = input(f'{PROMPTS[2]}\n')
        while not os.path.isfile(ref_genome_path) or not any(ref_genome_path.endswith(extension) for extension in ALLOWED_REF_GENOME_EXTENSIONS):
            ref_genome_path = input(f'\n\033[31m{PROMPTS[5]}\033[0m\n')

        print(f'\n\033[0m{PROMPTS[3]}\033[0m')
        ref_genome_tag = input(f'\n\033[31m{PROMPTS[4]}\033[0m\n\n').strip()
        while not ref_genome_tag:
            ref_genome_tag = input(f'\n\033[31m{PROMPTS[6]}\033[0m\n\n').strip()

        reference_genomes[ref_genome_path] = ref_genome_tag

        with open(ref_genomes_py, 'w') as f:
            f.write(f'reference_genomes = {reference_genomes}')

        ask_2_add_another_genome = input('Do you want to add another reference genome? (Y or N)\n').strip().lower()
        if ask_2_add_another_genome == 'n':
            break
        elif ask_2_add_another_genome != 'y':
            print('Invalid input!!')
            ask_2_add_another_genome = input('Do you want to add another reference genome? (Y or N)\n').strip().lower()

        else:
            continue

pick_ref_genome(add_ref_genome=True)
