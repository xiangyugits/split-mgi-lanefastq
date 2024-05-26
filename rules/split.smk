localrules:get_lane_index,lane_qc_stat


def Drc(sequence):
    comp_dict = {
        "A":"T",
        "T":"A",
        "G":"C",
        "C":"G",
        "a":"t",
        "t":"a",
        "g":"c",
        "c":"g",
        "N":"N",
    }
    #求互补序列
    sequence_list = list(sequence[::-1])
    sequence_list = [comp_dict[base] for base in sequence_list]
    string = ''.join(sequence_list)
    return string

def Drc2(sequence):
    if len(sequence)==16:
        return Drc(sequence[:8])+Drc(sequence[8:])
    elif len(sequence)==8:
        return Drc(sequence)
    else:
        return sequence

def get_lane_rc_index(info,lane,outfile,dualIndex=True):
    df=info[info.Lane==lane][["sampleName",'i5','i7']]
    
    if dualIndex==True:
        data=df.apply(lambda x:"".join([x.i5.strip(),x.i7.strip()]),axis=1).to_frame("Index")
    
    else:
        data=df.apply(lambda x:x.i7.strip(),axis=1).to_frame("Index")

    data.loc[:,'Name']=df.sampleName.str.replace(" ",'-')
    assert not data.Index.duplicated().any(), "包含重复序列"
    data.Index=data.Index.apply(Drc2)
    assert not (data.Index=="").any(), "样品index异常"
    data[["Name",'Index']].to_csv(outfile,index=False,header=False,sep='\t')

    return 0



rule get_lane_index:
    input:
        config["index_file"]
    output:
        join(temp_dir,"{Lane}.index_rc.csv")
    params:
        config["cycletype"]
    run:
        rawinfo=pd.read_csv(input[0],dtype={"sampleName":str})
        lane=os.path.basename(output[0]).split(".")[0]
        if config["cycletype"].count("r")==2:
            dualIndex=True
        else:
            dualIndex=False
        get_lane_rc_index(rawinfo,lane,output[0],dualIndex)#dualIndex false 时，仅采用i7


def get_split_input(wildcards):
    if is_paired_end:
        return dict(
            zip(
                ['fq1','fq2'],
                expand(
                    join(input_dir,"{Lane}",config["ChipID"]+"_{Lane}_read_{read}.fq.gz"),read=["1",'2'],**wildcards
                )
            )
        )
    else:
        return {'fq1':join(input_dir,"{Lane}/{chip}_{Lane}_read.fq.gz".format(chip=config["ChipID"],Lane=wildcards))}



rule run_splitBarcode:
    input:
        unpack(get_split_input),
        lane_index=join(temp_dir,"{Lane}.index_rc.csv")
    output:
        BarcodeStat=join(temp_dir,"{Lane}",'BarcodeStat.txt'),
    log:
        "log/{Lane}.log"
    threads:6
    cache:True
    params:
        cycle=config["cycletype"]
    script:
        "splitBarcode.py"

def get_stat(barcodestat,lane_summary):
    chip=barcodestat.split("/")[-3]
    lane=barcodestat.split("/")[-2]
    barcodestat=pd.read_csv(barcodestat,sep='\t')
    tmp1=pd.read_csv(lane_summary,index_col=0).iloc[[3,4,5,6,7],0]
    tmp2=pd.Series([chip,lane,barcodestat.iloc[-1,-2],barcodestat.iloc[-1,-1]],index=['Chip','Lane','TotalReads','SplitRate'])
    stat=pd.concat([tmp2,tmp1])
    stat=stat[['Chip','Lane','ChipProductivity(%)','CycleNumber','ImageArea','Q30(%)','SplitRate','TotalReads']]
    return stat

def get_BarcodeStats():
    if config["ChipID"][0]=="V":
        lane_list=["L01",'L02','L03','L04']

    else:
        lane_list=["L01","L02"]

    lane_range=config["lane_range"].split(",")

    lane_list= list(set(lane_list) & set(lane_range))


    split_output = expand(
        join(temp_dir,"{Lane}",'BarcodeStat.txt'),Lane=lane_list    
        )
    return split_output

def get_split_output():
    outfiles=[]
    for fq in fq_info.temp_fq.values:
        outfiles.append(fq)
        fqstat=fq.replace("gz",'fqStat.txt')
        outfiles.append(fqstat)
    return outfiles

rule lane_qc_stat:
    input:
        get_BarcodeStats()
    output:
        join(temp_dir,"{chip}.stat.txt".format(chip=config["ChipID"])),
    run:
        stat_out=[]
        for barcodestat in input:
            lane=barcodestat.split("/")[-2]
            lane_summary=join(input_dir,lane,'summaryTable.csv')
            tmp=get_stat(barcodestat,lane_summary)
            stat_out.append(tmp)
        df=pd.concat(stat_out,axis=1).T
        df.Chip=config["ChipID"]
        df.SplitRate=df.SplitRate.apply(lambda x:"{:.2f}".format(x))
        df.TotalReads=df.TotalReads.apply(lambda x:"{:.2f}".format(x/1e6))
        df=df.rename(columns={"SplitRate":'SplitRate(%)',"TotalReads":"TotalReads(M)"})
        df.sort_values(["Lane"]).to_csv(output[0],sep='\t',index=False)





