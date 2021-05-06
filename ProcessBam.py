import pysam
import subprocess as sp


# list of bam files
bamFiles = ["91KLH11.q30.autosome.bam", "91KLH18.q30.autosome.bam",
            "91KLM2.q30.autosome.bam", "BLSM27S.q30.autosome.bam",
            "BLSM41.q30.autosome.bam", "BLSM45.q30.autosome.bam",
            "DCZM17IV.q30.autosome.bam", "DCZ-M21II.q30.autosome.bam",
            "DCZ-M22IV.q30.autosome.bam", "DCZ-M6.q30.autosome.bam",
            "EDM124.q30.autosome.bam", "EDM139.q30.autosome.bam",
            "EDM176.q30.autosome.bam", "HJTM107.q30.autosome.bam",
            "HJTM109.q30.autosome.bam", "HJTM115.q30.autosome.bam",
            "HJTW13.q30.autosome.bam", "HMF32.q30.autosome.bam",
            "JCKM1-1.q30.autosome.bam", "JXNTM23.q30.autosome.bam",
            "JXNTM2.q30.autosome.bam", "LGM41.q30.autosome.bam",
            "LGM79.q30.autosome.bam", "LJM14.q30.autosome.bam",
            "LJM25.q30.autosome.bam", "LJM2.q30.autosome.bam",
            "LJM3.q30.autosome.bam", "LJM4.q30.autosome.bam",
            "LJM5.q30.autosome.bam", "MGS-M6.q30.autosome.bam",
            "MGS-M7L.q30.autosome.bam", "MGS-M7R.q30.autosome.bam",
            "MZGM10-1.q30.autosome.bam", "MZGM16.q30.autosome.bam",
            "MZGM25-2.q30.autosome.bam", "PLTM310.q30.autosome.bam",
            "PLTM311.q30.autosome.bam", "PLTM312.q30.autosome.bam",
            "PLTM313.q30.autosome.bam", "SM-SGDLM27.q30.autosome.bam",
            "SM-SGDLM6.q30.autosome.bam", "SM-SGDLM7X.q30.autosome.bam",
            "WD-WT1H16.q30.autosome.bam", "WD-WT5M2.q30.autosome.bam",
            "WGH35-1.q30.autosome.bam", "WGM20.q30.autosome.bam",
            "WGM35.q30.autosome.bam", "WGM43.q30.autosome.bam",
            "WGM70.q30.autosome.bam", "WGM76S.q30.autosome.bam",
            "WGM94.q30.autosome.bam", "WQM4.q30.autosome.bam",
            "XW-M1R18.q30.autosome.bam", "ZLNR-1.q30.autosome.bam",
            "ZLNR-2.q30.autosome.bam"]

for file in bamFiles:
    length = len(file)
    fileName = f"/home/classes/myanglab/data/earlyCN/{file}"
    sortedFileName = f"/home/classes/myanglab/data/earlyCN/{file[:length - 4]}.sorted.bam"
    sp.Popen(["samtools", "sort", fileName, "-o", sortedFileName], stdout = sp.PIPE)
    sp.Popen(["samtools", "index", sortedFileName], stdout = sp.PIPE)
    #pysam.sort("-o", file, sortedFileName) # length - 5 trims off the .bam
    #pysam.index(sortedFileName) I hate you pysam. Why won't you work?