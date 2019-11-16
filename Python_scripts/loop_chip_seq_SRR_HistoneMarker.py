import os

files = ["SRR577430_SRR577429.fastq","SRR227505_SRR227506.fastq","SRR577392_SRR577393.fastq","SRR227473_SRR227472.fastq","SRR227445_SRR227446.fastq","SRR568261_SRR568260.fastq","SRR350914_SRR350915.fastq","SRR227556_SRR227557.fastq","SRR5338600_SRR5338596_SRR5338597_SRR5338598_SRR5338599.fastq","SRR577378_SRR577379.fastq","SRR227441_SRR227442.fastq","SRR227456_SRR227457.fastq","SRR568255_SRR568254.fastq","SRR227602_SRR227603.fastq","SRR227413_SRR227414.fastq"]
labels = ["H3K36me3","H3K36me3","H3K27me3","H3K27me3","H3K27ac","H2AFZ","H3K4me1","H3K4me2","H3K4me3","H3K4me3","H3K4me3","H4K20me1","H3K9me3","H3K9ac","H3K79me2"]

for k in range(0,len(files)):
    #os.system("bash /cta/users/ardacetin/globalRepair/chipseq/partial_chipseq.sh " + files[k] + " " + labels[k])
    print("bash /cta/users/ardacetin/globalRepair/chipseq/partial_chipseq.sh " + files[k] + " " + labels[k])
