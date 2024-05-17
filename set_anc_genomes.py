import os
from VCSeek_config import PROMPTS, ALLOWED_REF_GENOME_EXTENSIONS

REF_GENOMES = {}

ref_genome_path = input(f'{PROMPTS[2]}\n')
while not os.path.isfile(ref_genome_path) and any(ref_genome_path.endswith(extension) for extension in ALLOWED_REF_GENOME_EXTENSIONS):
    ref_genome_path = input(f'\n\033[31m{PROMPTS[5]}\033[0m\n')

print(f'\n\033[0m{PROMPTS[3]}\033[0m')
ref_genome_tag = input(f'\n\033[31m{PROMPTS[4]}\033[0m\n\n').strip()
while not ref_genome_tag:
    ref_genome_tag = input(f'\n\033[31m{PROMPTS[6]}\033[0m\n\n').strip()

REF_GENOMES[ref_genome_path] = ref_genome_tag
print(REF_GENOMES)
