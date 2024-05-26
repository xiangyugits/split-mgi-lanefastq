__author__ = "zhouxiangyu"
__copyright__ = ""
__email__ = "zhouxiangyu@cibr.ac.cn"
__license__ = "CIBR"


import os
#import tempfile
from snakemake.shell import shell


extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)


fq1 = snakemake.input.get("fq1")
assert fq1 is not None, "input-> fq1 is a required input parameter"


fq2 = snakemake.input.get("fq2")
if fq2:
    assert len(fq1) == len(fq2), "input-> equal number of files required for fq1 and fq2"

input_str_fq1 = "-1 "+fq1
input_str_fq2 = "-2 "+fq2 if fq2 is not None else ""
input_str = " ".join([input_str_fq1, input_str_fq2])


cycle = snakemake.params.get("cycle", "r100e1r100e1b8b8")

cycle2paramdict={
    "r100e1r100e1b8b8":"-b 200 8 1 -b 208 8 1",
    "r100e1r100e1b8":"-b 208 8 1",
    "r150e1r150e1b8b8":"-b 300 8 1 -b 308 8 1",
    "r50e1b8":"-b 50 8 1"   
}

lane_index=snakemake.input.get("lane_index")
cycleparam=cycle2paramdict[cycle]
outdir=os.path.split(snakemake.output.BarcodeStat)[0]

#with tempfile.TemporaryDirectory() as tmpdir:
shell(
    "ulimit -n 10000;"
    "/home/zhangli_lab/zhouxiangyu/DATA/software/splitBarcode_lite/bin/splitBarcode"
    " -B {lane_index}"
    " {input_str}" 
    " -o {outdir}"
    " {cycleparam}"
    " -t {snakemake.threads}"
)