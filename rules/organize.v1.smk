
localrules:get_SAMPLE_fq,fastq_stat,sample_stat

fq_info=get_fq_info(config["index_file"],is_paired_end,output_dir,temp_dir)


# def get_sample_stats(wildcards):
#     filter_info=fq_info.loc[(fq_info.lab_user==wildcards.Lab_user),]
#     return filter_info.temp_fq.str.replace("gz","fqStat.txt").to_list()

def get_sample_stats(wildcards):
    stat = checkpoints.lane_qc_stat.get(**wildcards).output[0]
    fq_output=os.path.split(stat)[0]
    SMP,= glob_wildcards(join(fq_output,"{sample}.fq.fqStat.txt"))
    return expand(
        join(fq_output,"{SMP}.fq.fqStat.txt"),SMP=SMP
    )


rule fastq_stat:
    input:
        get_sample_stats
    output:
        join(output_dir,'{Lab_user}',"fq.QC.stat.csv")

    run:
        res_list=[]
        for i in input:
            tmp=read_fq_stat(i)
            res_list.append(tmp)
        res=pd.concat(res_list)

        res=fq_info.merge(res,left_on="temp_fq_name",right_index=True)
        res[["ReadNum",'row_readLen','BaseNum']]=res[["ReadNum",'row_readLen','BaseNum']].astype(int)
        res[["GC%",'Q10%','Q20%','Q30%']]=res[["GC%",'Q10%','Q20%','Q30%']].astype(float)
        res1=res.groupby("final_fq")[["ReadNum",'BaseNum']].sum()
        res2=res.groupby("final_fq")[["row_readLen",'GC%','Q10%','Q20%','Q30%']].mean().astype(float)
        out=pd.concat([res1,res2],axis=1)

        out=fq_info[["final_fq"]].merge(out,left_on="final_fq",right_index=True).rename(columns={"final_fq":'fq','row_readLen':'ReadLength'})
        out.fq=out.fq.apply(lambda x:os.path.basename(x))
        out.ReadLength=out.ReadLength.astype(int)
        out=out.sort_index()
        out.to_csv(output[0])


rule sample_stat:
    input:
        join(output_dir,'{Lab_user}',"fq.QC.stat.csv")
    output:
        join(output_dir,'{Lab_user}',"Sample.QC.stat.txt")
    run:
        res=pd.read_csv(input[0])
        res=res.groupby("sampleName")[["ReadNum",'BaseNum',"ReadLength",'GC%','Q10%','Q20%','Q30%']].mean()
        out=fq_info[['batch']].merge(res,left_index=True,right_index=True,how='outer').fillna(0)
        out[["ReadNum",'ReadLength','BaseNum']]=out[["ReadNum",'ReadLength','BaseNum']].astype(int)
        out[["GC%",'Q10%','Q20%','Q30%']]=out[["GC%",'Q10%','Q20%','Q30%']].map(lambda x:"{:.2f}".format(x))
        out.reset_index().sort_values(["batch",'sampleName']).to_csv(output[0],sep='\t',index=False)


#######################################################################################################################

def get_SAMPLE_input(wildcards):
    stat = checkpoints.lane_qc_stat.get(**wildcards).output[0]
    fq_output=os.path.split(stat)[0]
    lanelist,temp_fq_list= glob_wildcards(join(fq_output,"{Lane}","{temp_fq_name}.fq.gz"))
    final_fq_info=fq_info.loc[fq_info.temp_fq_name.isin(temp_fq_list),]
    if is_paired_end:
        fq_list= final_fq_info.loc[(final_fq_info.SAMPLE==wildcards.SAMPLE)&(final_fq_info.read==wildcards.read),].temp_fq.to_list()
    else:
        fq_list= final_fq_info.loc[(final_fq_info.SAMPLE==wildcards.SAMPLE),].temp_fq.to_list()

    assert len(fq_list)>0
    return fq_list


def get_SAMPLE_output():
    if is_paired_end:
        return join(output_dir,'{SAMPLE}.{read}.fq.gz')
    else:
        return join(output_dir,'{SAMPLE}.fq.gz')


rule get_SAMPLE_fq:
    input:
        get_SAMPLE_input
    output:
        get_SAMPLE_output()        
    shell:
        "cat {input} > {output}"


def Group_fq_done(wildcards):
    stat = checkpoints.lane_qc_stat.get(**wildcards).output[0]
    fq_output=os.path.split(stat)[0]
    lanelist,temp_fq_list= glob_wildcards(join(fq_output,"{Lane}","{temp_fq_name}.fq.gz"))
    final_fq_info=fq_info.loc[fq_info.temp_fq_name.isin(temp_fq_list),]
    Group_info=final_fq_info[final_fq_info.lab_user==wildcards.Lab_user]
    #import ipdb;ipdb.set_trace()
    if is_paired_end:
        return expand(
            join(output_dir,'{SAMPLE}.{read}.fq.gz'),SAMPLE=Group_info.SAMPLE.unique(),read=["R1",'R2']
        )
    else:
        return expand(
            join(output_dir,'{SAMPLE}.fq.gz'),SAMPLE=Group_info.SAMPLE.unique()
       )


rule MD5:
    input:
        Group_fq_done
    output:
        join(output_dir,'{Lab_user}',"MD5.txt")
    shell:
        "dirpath=`dirname {output}`;"
        "cd $dirpath;"
        "md5sum `ls */*fq.gz` > MD5.txt"

#######################################################################################################################

def final_stat_output():
    return expand(
        join(output_dir,'{Lab_user}',"{stat_file}"),
        Lab_user=fq_info.lab_user.unique(),stat_file=["Sample.QC.stat.txt","MD5.txt"]
    )