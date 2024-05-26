import pandas as pd
import numpy as np
import sys,glob,os,warnings


def get_standard_index(csv):
    CHIP=csv.split(".")[0]
    excel=pd.read_csv(csv,index_col=0)
    excel=excel[excel.sampleName.isna()==False] #去除空值
    #重新整理lab/from 信息
    excel.belonguser=excel.apply(lambda x:x["name"] if (pd.isna(x.belonguser))|(x.belonguser=="unkown") else x.belonguser,axis=1)
    excel.belonggroup=excel.apply(lambda x:x.lab if (pd.isna(x.belonggroup))|(x.belonggroup=="unkown") else x.belonggroup,axis=1)
    excel=excel[['sampleName','i5','i7',"Lane","lab","belonguser","batch","libtype"]].rename(columns={'lab':"Lab"}).sort_values(['Lane','Lab','sampleName'])
    
    #i5/i7
    excel[['i5','i7']]=excel[['i5','i7']].fillna("NNNNNNNN")
    excel[['i5','i7']]=excel[['i5','i7']].applymap(lambda x:x.strip(" "))
    excel.i5=excel.i5.apply(lambda x:x[2:] if len(x)==10 else x)
    excel.i7=excel.i7.apply(lambda x:x[2:] if len(x)==10 else x)   
    
    #检查index重复
    excel.groupby("Lane").apply(check_dup)
    excel.libtype=excel.libtype.fillna("unkown")
    excel.insert(3,"CHIP",CHIP)
    return excel
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
    return Drc(sequence[:8])+Drc(sequence[8:])

def check_dup(df):
    df2=df[["Lane","i5",'i7']][df.duplicated()]
    if len(df2)>0:
        print(df2.Lane.values[0])
        print(df2)    
        return 1
    else:
        return 0

def info_zh_to_en(zh):
    info_dict={
        "罗敏敏实验室":"LuoLab","白凌实验室":"BaiLab",
        "陈坚实验室":"ChenLab","韩东实验室":"HanLab",
        "戈鹉平实验室":"GeLab","孙文智实验室":"SunLab",
        "梅林实验室":"MeiLab","王同飞实验室":"WangLab",
        "井淼实验室":"JingLab",
        
        "戴睿成":"DaiRC","付佳颖":"FuJY",
        "汪意":"WangYi","王雅":"WangYa","王睿宇":"WangRY",
        
    }
    if zh in info_dict.keys():
        return info_dict[zh]
    else:
        return zh

def smartseq_index(excel):
    excel.libtype="Smartseq"
    excel.sampleName=excel.apply(lambda x:x.batch+"-"+x.sampleName if "-" not in x.sampleName else x.sampleName,axis=1)
    excel.i7=excel.i7.apply(Drc)
    return excel

def split_info(Laneinfo,subnum,maxnum):
    Laneinfo=Laneinfo.reset_index(drop=True)
    if Laneinfo.shape[0]<maxnum:
        return Laneinfo
    Laneinfo.Lane=Laneinfo.apply(lambda x:x.Lane+"_"+str(int(x.name/subnum)+1),axis=1)
    return Laneinfo


def check_lane(lane,SampleInfo):
    df=pd.read_csv("temp_dir/%s/SequenceStat.txt"%lane,sep='\t')
    df.columns=pd.Series(df.columns).apply(lambda x:x.strip()).values
    df.Barcode=df.Barcode.apply(lambda x:x.strip())
    df=df[(df.Barcode=="undecoded")&(df.Count>1e5)]

    df.loc[:,'I5']=df.iloc[:,0].apply(lambda x:x[:8])
    df.loc[:,'I5_Drc']=df.iloc[:,0].apply(lambda x:Drc(x[:8]))
    df.loc[:,'I7']=df.iloc[:,0].apply(lambda x:x[8:])
    df.loc[:,'I7_Drc']=df.iloc[:,0].apply(lambda x:Drc(x[8:]))

    for i5 in ["I5",'I5_Drc']:
        for i7 in ["I7",'I7_Drc']:
            out=df.rename(columns={i5:"i5",i7:"i7"}).merge(SampleInfo,on=["i5",'i7'])
            out=out[["Lane","belonguser","batch"]].value_counts().sort_index()
            if out.shape[0]>0:
                print(i5,i7)
                print(out)
                
def get_stat(dirname):
    reslist=[]
    for i in glob.glob("{0}/*/Sample.QC.stat.txt".format(dirname)):
        tmp=pd.read_csv(i,sep='\t')
        user=i.split("/")[-2]
        tmp.loc[:,'user']=user
        reslist.append(tmp)
    res=pd.concat(reslist)
    return res
        
def get_stat_result(res,groupids):
    tmp=res.drop_duplicates()
    tmp1=tmp.groupby(groupids).ReadNum.sum()
    tmp2=tmp.groupby(groupids).sampleName.count()
    tmp3=tmp.groupby(groupids).BaseNum.sum()
    out=pd.concat([tmp1,tmp2,tmp3],axis=1,keys=['reads','Samples','Size'])
    out.loc[:,'Size(G)']=out.Size.apply(lambda x:"{:.2f}".format(x/1e9))
    out.loc[:,'reads']=out.reads.apply(lambda x:"{:.2f}".format(x/1e6))
    return out.sort_index(axis=1).rename(columns={"reads":"reads(M)"}).drop("Size",axis=1)