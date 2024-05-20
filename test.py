import os
import ast

dict_py = os.path.expanduser("/home/rmeneses/EVB_Variant_Calling/ref_genomes_db.py")

with open(dict_py, 'r') as f:
    content = f.read()
    tree = ast.parse(content, mode='exec')
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == "reference_genomes":
                    ast.literal_eval(node.value)