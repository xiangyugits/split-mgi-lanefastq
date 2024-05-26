import glob
from os.path import join
import pandas as pd

def read_fq_stat(file):
    with open(file,'r') as f:
        stat_info=f.readlines()[1:11]
        se=pd.Series(stat_info).str.strip().str.strip("#").iloc[[1,2,4,6,7,8,9]].str.split("\t",expand=True).set_index(0).iloc[:,0]
        df=se.to_frame(file.split("/")[-1].split(".fq.fqStat.txt")[0])
        return df.T

def get_fq_info(index_file,is_paired_end,outdir,tempdir):
    info=pd.read_csv(index_file,dtype={"sampleName": str}).set_index("sampleName")
    info.loc[:,'outdir']=outdir
    info.loc[:,'tempdir']=tempdir
    info.loc[:,"SAMPLE"]=info.apply(lambda x:"/".join([x.Lab,x.belonguser,x.batch,x.name]),axis=1)
    info.loc[:,'fq_dir']=info.apply(lambda x:join(x.outdir,x.Lab,x.belonguser,x.batch),axis=1)
    if is_paired_end:
        info1=info.copy()
        info1.loc[:,'temp_fq']= info1.apply(lambda x:join(x.tempdir,x.Lane,"{0}_{1}_{2}_1.fq.gz".format(x.CHIP,x.Lane,x.name)),axis=1)
        info1.loc[:,'final_fq']=info1.apply(lambda x:join(x.fq_dir,x.name+"_R1.fq.gz"),axis=1)
        info1.loc[:,'read']="R1"
        info2=info.copy()
        info2.loc[:,'temp_fq']=info2.apply(lambda x:join(x.tempdir,x.Lane,"{0}_{1}_{2}_2.fq.gz".format(x.CHIP,x.Lane,x.name)),axis=1)
        info2.loc[:,'final_fq']=info2.apply(lambda x:join(x.fq_dir,x.name+"_R2.fq.gz"),axis=1)
        info2.loc[:,'read']="R2"
        fq_info=pd.concat([info1,info2])
    else:
        info1=info.copy()
        info1.loc[:,'temp_fq']=info1.apply(lambda x:join(x.tempdir,x.Lane,"{0}_{1}_{2}.fq.gz".format(x.CHIP,x.Lane,x.name)),axis=1)
        info1.loc[:,'final_fq']=info1.apply(lambda x:join(x.fq_dir,x.name+".fq.gz"),axis=1)
        info1.loc[:,'read']="read"
        fq_info=info1
    
    fq_info.loc[:,'lab_user']=fq_info.apply(lambda x:"/".join([x.Lab,x.belonguser]),axis=1)
    fq_info.loc[:,'Lab_user_batch']=fq_info.apply(lambda x:"/".join([x.Lab,x.belonguser,x.batch]),axis=1)
    fq_info.loc[:,'temp_fq_name']=fq_info.apply(lambda x:os.path.basename(x.temp_fq).split(".fq.gz")[0],axis=1)
    fq_info.loc[:,'sample_read']=fq_info.apply(lambda x:os.path.basename(x.final_fq).split(".fq.gz")[0],axis=1)

    fq_info=fq_info[['SAMPLE',"lab_user","Lab_user_batch",'sample_read','temp_fq_name',"temp_fq",'final_fq','read','batch']]
    return fq_info
