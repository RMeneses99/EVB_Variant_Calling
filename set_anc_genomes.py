import os
import ast
from VCSeek_config import PROMPTS, ALLOWED_REF_GENOME_EXTENSIONS

def pick_ancestral_genome(add_ancestral=False):
    script_path = os.path.dirname(os.path.realpath(__file__))
    # Path to the file containing reference genomes
    anc_file = os.path.join(script_path,'ancestral_genomes_db.py')
    # Check if the reference file exists
    ref_exists = os.path.isfile(anc_file)
    ancestral_genomes = {}
    # Load existing reference genomes if the file exists
    if ref_exists:
        with open(anc_file, 'r') as file:
            content = file.read()
            tree = ast.parse(content, mode='exec')
            for node in tree.body:
                if isinstance(node, ast.Assign):
                    for target in node.targets:
                        if isinstance(target, ast.Name) and target.id == "ancestral_genomes":
                            ancestral_genomes = ast.literal_eval(node.value)

        while add_ancestral:
            # Get path for new reference genome
            ancestral_path = input('Enter the path to the reference genome:\n')
            while not os.path.isfile(ancestral_path) or not any(ancestral_path.endswith(ext) for ext in ALLOWED_REF_GENOME_EXTENSIONS):
                ancestral_path = input('\nInvalid path or extension. Enter a valid path to the ancestral reference genome:\n')

            # Get tag for the new reference genome
            ancestral_tag = input('Enter a tag for the ancestral reference genome:\n').strip()
            while not ancestral_tag:
                ancestral_tag = input('Tag cannot be empty. Enter a tag for the reference genome:\n').strip()
        
            ancestral_genomes[ancestral_path] = ancestral_tag 

            # Ask if the user wants to add another reference genome
            add_another = input('Do you want to add another ancestral reference genome? (Y or N)\n').strip().lower()
            while add_another not in ['y','n']:
                print('Invalid input!!\n')
                add_another = input('Do you want to add another reference genome? (Y or N)\n').strip().lower()
            if add_another == 'n':
                add_ancestral == False

            # Write the updated reference genomes to the file
            with open(anc_file, 'w') as file:
                file.write(f'ancestral_genomes = {ancestral_genomes}')
            add_ancestral = False

    else:
         # Get path for new reference genome
        ancestral_path = input('Enter the path to the reference genome:\n')
        while not os.path.isfile(ancestral_path) or not any(ancestral_path.endswith(ext) for ext in ALLOWED_REF_GENOME_EXTENSIONS):
            ancestral_path = input('\nInvalid path or extension. Enter a valid path to the ancestral reference genome:\n')

        # Get tag for the new reference genome
        ancestral_tag = input('Enter a tag for the ancestral reference genome:\n').strip()
        while not ancestral_tag:
            ancestral_tag = input('Tag cannot be empty. Enter a tag for the reference genome:\n').strip()
        
        ancestral_genomes[ancestral_path] = ancestral_tag

        # Ask if the user wants to add another reference genome
        add_another = input('Do you want to add another ancestral reference genome? (Y or N)\n').strip().lower()
        while add_another == 'y':
                # Get path for new reference genome
            ancestral_path = input('Enter the path to the reference genome:\n')
            while not os.path.isfile(ancestral_path) or not any(ancestral_path.endswith(ext) for ext in ALLOWED_REF_GENOME_EXTENSIONS):
                ancestral_path = input('\nInvalid path or extension. Enter a valid path to the ancestral reference genome:\n')

            # Get tag for the new reference genome
            ancestral_tag = input('Enter a tag for the ancestral reference genome:\n').strip()
            while not ancestral_tag:
                ancestral_tag = input('Tag cannot be empty. Enter a tag for the reference genome:\n').strip()
        
            ancestral_genomes[ancestral_path] = ancestral_tag
            # Ask if the user wants to add another reference genome
            add_another = input('Do you want to add another ancestral reference genome? (Y or N)\n').strip().lower()
        
        while add_another not in ['y','n']:
            print('Invalid input!!\n')
            add_another = input('Do you want to add another reference genome? (Y or N)\n').strip().lower()

        if add_another == 'n':
            with open(anc_file, 'x') as file:
                file.write(f'ancestral_genomes = {ancestral_genomes}')



def main():
    pick_ancestral_genome()

if __name__ == '__main___':
    main()