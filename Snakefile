
from os.path import join
import glob

#### init ####
input_dir=join(config["rawfq_dir"],config["ChipID"])
temp_dir = config["temp_dir"]
output_dir= config["outfq_dir"] 

include: "rules/common.smk"


is_paired_end= (config["cycletype"].count("r")==2)
#samples=pd.read_csv(config['index_file'],dtype={"sampleName": str}).set_index("sampleName")
fq_info=get_fq_info(config["index_file"],is_paired_end,output_dir,temp_dir)
Step=config["step"]
#import ipdb;ipdb.set_trace()


include: "rules/split.smk"
include: "rules/organize.smk"


def final_output(Step):
    
    if float(Step)==1:
        return join(temp_dir,"{chip}.stat.txt".format(chip=config["ChipID"]))
    else:
        return final_stat_output()
rule all:
    input:
        final_output(Step)






        







    