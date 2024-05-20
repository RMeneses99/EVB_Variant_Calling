import Toolkit

path_2_analyze = Toolkit.ask_usr(breseq=True)
#Toolkit.quali_check(path_2_analyze)
read_pairs = Toolkit.pair_by(path_2_analyze)
Toolkit.pre_processing(path_2_analyze, hgt=True)




