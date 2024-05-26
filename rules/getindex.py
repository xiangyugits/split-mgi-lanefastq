#!/home/zhangli_lab/zhouxiangyu/DATA/anaconda3/bin/python
#-*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import os,time,glob,argparse,sys
from sqlalchemy import create_engine

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="""

""")
parser.add_argument('-i', dest='input', help='chipID', required=True)

args = parser.parse_args()
ID=args.input


def main(chip):
    engine=create_engine("mysql+pymysql://root:cibr123!@117.133.167.134:8138/infoseq")
    chipinfo =pd.read_sql_query("select * from cibr_sys_task WHERE chipno='%s'"%chip, engine)
    chipid=chipinfo.id.values[0]
    lanenum=4 if ("V" in chip)|("v" in chip) else 2
    laneids=[chipid+str(num) for num in range(1,lanenum+1)]
    laneinfo = pd.DataFrame()
    for i in laneids:
        tmp=pd.read_sql_query("select * from cibr_sys_info WHERE taskId='%s'"%i, engine)
        laneinfo=laneinfo.append(tmp)

    laneinfo= laneinfo[['id','name','lab','batch']].rename(columns={"id":"infoId"})

    data=pd.read_sql_query("select * from cibr_sys_index WHERE taskId='%s'"%chipid, engine)
    data2=data[['sampleName','i5','i7','libtype','belonguser','belonggroup','infoId']].copy()
    data2.loc[:,'Lane']=data.laneId.apply(lambda x:"L0"+x[-1])
    result=data2.merge(laneinfo,on='infoId').drop("infoId",axis=1)
    result.to_csv("{}.index.csv".format(chip))


if __name__=='__main__':
    main(ID)
    

