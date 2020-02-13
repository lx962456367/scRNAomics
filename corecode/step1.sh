# 小鼠单细胞多组学分析流程

=====================================================================================================================
第一步：先搞到数据（通过已发表文献/数据库下载成对的RNA测序数据：scRNA和bulkRNA）

# bulk-RNA TPM result（晓瑜分析bulkRNA数据）
#肝
/dellfsqd2/ST_LBI/USER/zhouxiaoyu/SeqTech/liver/output/TPM_liver.txt
#肾
/dellfsqd2/ST_LBI/USER/zhouxiaoyu/SeqTech/kidney/output/kidney_TPM.txt
#肺
/dellfsqd2/ST_LBI/USER/zhouxiaoyu/SeqTech/lung/output/lung_TPM.txt
#心
/dellfsqd2/ST_LBI/USER/zhouxiaoyu/SeqTech/heart/output/heart_TPM.tx
#乳腺
/dellfsqd2/ST_LBI/USER/zhouxiaoyu/SeqTech/breast/output/breast_TPM.txt
#脾
/dellfsqd2/ST_LBI/USER/zhouxiaoyu/SeqTech/spleen/output/spleen_TPM.txt
#胸腺
/dellfsqd2/ST_LBI/USER/zhouxiaoyu/SeqTech/thymus/output/thymus_TPM.txt
#结肠
/dellfsqd2/ST_LBI/USER/zhouxiaoyu/SeqTech/colon/output/colon_TPM.txt



# single-cell RNA raw_matrix (barcodes.tsv  genes.tsv  matrix.mtx) 
#肝
/dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/droplet/Liver-10X_P4_2
#肾
/dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/droplet/Kidney-10X_P7_5
#肺
/dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/droplet/Lung-10X_P7_9
#心
/dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/droplet/Heart-10X_P7_4
#乳腺
/dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/droplet/Mammary-10X_P7_12
#脾
/dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/droplet/Spleen-10X_P4_7
#胸腺
/dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/droplet/Thymus-10X_P7_11

$ ls /dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/droplet/Liver-10X_P4_2/
barcodes.tsv  genes.tsv  matrix.mtx

#结肠
/dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/FACS/Colon-counts.csv

$ less -S Colon-counts.csv
# "","F11.MAA001871.3_39_F.1.1","G16.MAA001871.3_39_F.1.1","I1.MAA001871.3_39_F.1.1","A13.MAA001871.3_39_F.1.1","B16.MAA001871.3_39_F.1.1","C18.MAA001871.3_39_F
# "0610005C13Rik",0,0,0,0,0,0,0,0,0,0,3,80,22,0,23,12,21,57,0,0,0,0,0,33,0,0,0,0,1,0,0,13,0,0,61,125,0,3,0,0,34,0,0,47,68,0,0,0,3,0,0,0,6,0,0,0,0,1,0,0,0,0,0,0,
# "0610007C21Rik",0,0,0,0,0,0,0,0,361,154,0,114,0,39,51,1,191,148,0,0,0,0,80,82,0,0,0,0,38,0,175,2,104,0,0,249,1,171,0,44,14,0,236,117,140,180,110,0,20,1,40,19,
# "0610007L01Rik",0,122,44,0,3,0,22,0,179,0,93,0,0,28,0,0,248,122,0,0,1,0,60,59,0,0,1,0,0,0,0,26,0,0,0,0,0,90,0,0,0,0,154,16,0,0,138,0,0,6,0,0,6,0,284,0,0,59,0,
# "0610007N19Rik",1,0,0,2,0,3,1,4,74,121,0,0,0,1,5,0,87,2,0,0,2,0,26,0,0,2,1,0,2,0,1,0,2,0,0,22,1,92,1,0,0,2,124,79,34,41,54,9,0,0,0,5,8,0,66,0,0,71,0,6,46,7,0,

# single-cell RNA artifical_matrix (barcodes.tsv  genes.tsv  matrix.mtx) 
#肝
/dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/artifical_matrix/Liver/Liver-10X_P4_2.csv 
#肾
/dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/artifical_matrix/Kidney/Kidney-10X_P7_5.csv 
#肺
/dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/artifical_matrix/Lung/Lung-10X_P7_9.csv 
#心
/dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/artifical_matrix/Heart/Heart-10X_P7_4.csv 
#乳腺
/dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/artifical_matrix/Mammary/Mammary-10X_P7_12.csv 
#脾
/dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/artifical_matrix/Spleen/Spleen-10X_P4_7.csv 
#胸腺
/dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/artifical_matrix/Thymus/Thymus-10X_P7_11.csv 


=====================================================================================================================
第二步：对单细胞数据进行预处理（将下载的原始稀疏矩阵-->容易理解的密集矩阵）

# 使用cellranger mat2csv命令将原始稀疏矩阵转换为密集CSV格式
$ cd /dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/artifical_matrix/
$ qsub -cwd -l vf=1G -l num_proc=3 -q st.q -P P18Z19700N0012 run_cellranger_mat2csv.sh
$ head run_cellranger_mat2csv.sh
#/share/biosoft/pipeline/Single_Cell/Software/cellranger-3.0.2_R/cellranger mat2csv /dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/droplet/Bladder-10X_P4_3 Bladder-10X_P4_3.csv 
#/share/biosoft/pipeline/Single_Cell/Software/cellranger-3.0.2_R/cellranger mat2csv /dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/droplet/Bladder-10X_P4_4 Bladder-10X_P4_4.csv 
#/share/biosoft/pipeline/Single_Cell/Software/cellranger-3.0.2_R/cellranger mat2csv /dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/droplet/Bladder-10X_P7_7 Bladder-10X_P7_7.csv 
#/share/biosoft/pipeline/Single_Cell/Software/cellranger-3.0.2_R/cellranger mat2csv /dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/droplet/Heart-10X_P7_4 Heart-10X_P7_4.csv 
#/share/biosoft/pipeline/Single_Cell/Software/cellranger-3.0.2_R/cellranger mat2csv /dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/droplet/Kidney-10X_P4_5 Kidney-10X_P4_5.csv 
#/share/biosoft/pipeline/Single_Cell/Software/cellranger-3.0.2_R/cellranger mat2csv /dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/droplet/Kidney-10X_P4_6 Kidney-10X_P4_6.csv 
#/share/biosoft/pipeline/Single_Cell/Software/cellranger-3.0.2_R/cellranger mat2csv /dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/droplet/Kidney-10X_P7_5 Kidney-10X_P7_5.csv 
#/share/biosoft/pipeline/Single_Cell/Software/cellranger-3.0.2_R/cellranger mat2csv /dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/droplet/Liver-10X_P4_2 Liver-10X_P4_2.csv 
#/share/biosoft/pipeline/Single_Cell/Software/cellranger-3.0.2_R/cellranger mat2csv /dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/droplet/Liver-10X_P7_0 Liver-10X_P7_0.csv 
#/share/biosoft/pipeline/Single_Cell/Software/cellranger-3.0.2_R/cellranger mat2csv /dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/droplet/Liver-10X_P7_1 Liver-10X_P7_1.csv 

# output
$ head Heart-10X_P7_4.csv
# Gm5506,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,2,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0
# 1700044K03Rik,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
# Dtwd2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,
# 1700065O20Rik,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
# Dmxl1,0,0,0,2,0,0,0,0,0,1,0,1,0,0,0,0,0,1,0,0,0,1,0,1,0,0,0,0,0,0,2,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
# Tnfaip8,0,0,0,1,0,0,0,1,0,0,0,0,0,0,3,3,0,0,1,0,2,1,3,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,1,0,0,0,1,0,3,0,0,0,1,0,0,0,0,0,2,0,0,0,0,2,0,0,0,
# ......

=====================================================================================================================
第三步：对CSV格式的矩阵进行横行求和，计算的是每个基因在所有细胞中的覆盖总数(以小鼠“乳腺”mammary为例)

# awk对行进行求和（已跳过首行的细胞名称，是从第二行开始计算的）
$ cd /dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/artifical_matrix/Mammary/
$ awk -F"," '{b[NR]=$0; for(i=3;i<=NF;i++)a[NR]+=$i;}END{for(i=2;i<=NR;i++) print a[i]}' Mammary-10X_P7_12.csv > Mammary-10X_P7_12_add.txt
# output
$ head Mammary-10X_P7_12_add.txt
# 1
# 0
# 304
# 831
# 421
# 1405
# 49
# 390
# 0
# 0

# 从矩阵中提取首列——基因名称
$ awk -F"," '{print$1}' Mammary-10X_P7_12.csv > Mammary-10X_P7_12_header.txt
#output
$ head Mammary-10X_P7_12_header.txt
#Xkr4
#Rp1
#Sox17
#Mrpl15
#Lypla1
#Tcea1
#Rgs20
#Atp6v1h
#Oprk1
#Npbwr1

# 将xx_add.txt和xx_header.txt合并到一个文件中
$ paste Mammary-10X_P7_12_header.txt Mammary-10X_P7_12_add.txt > Mammary-10X_P7_12_final.txt
#output
$ head Mammary-10X_P7_12_final.txt
Xkr4	1
Rp1	0
Sox17	304
Mrpl15	831
Lypla1	421
Tcea1	1405
Rgs20	49
Atp6v1h	390
Oprk1	0
Npbwr1	0

=====================================================================================================================
第四步：计算RPM值

# 先通过sftp将scRNA的输出文件get到本地指定文件夹中；
sftp> get /dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/artifical_matrix/Mammary/Mammary-10X_P7_12_final.txt ./Desktop

# 然后,用excel表打开Mammary-10X_P7_12_final.txt文件
# 新增四列，分别是：
SYMBOL, total reads, All reads (x), RPM (B value)
基因名称，每个基因的reads总数（所有细胞都算），第二列加和，第二列先乘以1000000再除以第三列

=====================================================================================================================
第五步：有了RPM值之后，接下来就TPM和RPM取交集

# 先通过sftp将bulk RNA的输出文件get到本地指定文件夹中；
sftp> get /dellfsqd2/ST_LBI/USER/zhouxiaoyu/SeqTech/spleen/output/mammary_TPM.txt ./Desktop

# 然后,用excel表打开mammary_TPM.txt文件，并用Average()函数计算出每个基因在不同样本之间TPM的均值；
# 在同一个excel表中新建一个sheet1，复制TPM均值那一列，选择性粘贴"值"到sheet1中；
# 在sheet1中，用vlookup()函数，从RPM所在的excel表中抓取与TPM有着相同基因名称的RPM值,并通过筛选功能去掉NA的行；（剩下的是交集基因）
# 在sheet2中，用log()函数，对TPM和RPM两列分别取log(X+1,2)；
# 在sheet3中，复制log2之后的这两列，选择性粘贴值到sheet3中；


第六步：做相关图
# 复制sheet3中这两列，粘贴到Prism模板（Collection_Mammray.pzfx）中做相关性散点图。

=====================================================================================================================
# Start anaysis 
$ mkdir -p /dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/seurat_liver 
$ source /dellfsqd2/ST_LBI/USER/panxiaoguang/app/miniconda3/bin/activate Rlibs
$ R
library(Seurat)
library(cowplot)
library("dplyr")
library(reticulate)
library(ggplot2)
use_condaenv("r-reticulate")

## Liver
# Load the liver dataset
liver <- Read10X(data.dir = '/dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/droplet/Liver-10X_P4_2')
liver <- CreateSeuratObject(counts = liver)
liver
# An object of class Seurat 
# 23433 features across 1006 samples within 1 assay 
# Active assay: RNA (23433 features)

#liver[["percent.mt"]] <- PercentageFeatureSet(object = liver, pattern = "^mt-")
#pdf("liver_QC.pdf");VlnPlot(liver, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3);dev.off()
# Normalize the data
liver <- NormalizeData(liver)

# Identify highly variable features
liver <- FindVariableFeatures(liver, selection.method = "vst", nfeatures = 2000)

# Scale the data
all.genes <- rownames(liver)
liver <- ScaleData(liver, features = all.genes)

# Perform PCA on the scaled data
#liver <- RunPCA(liver, features = FindVariableFeatures(object = liver))
liver <- RunPCA(liver, npcs = 30, verbose = FALSE)

# Cluster the cells
liver <- RunUMAP(liver, reduction = "pca", dims = 1:20)
liver <- FindNeighbors(liver, reduction = "pca", dims = 1:20)
liver <- FindClusters(liver, resolution = 1.0); table(liver@active.ident)
#  0   1   2   3   4   5   6   7   8 
# 232 193 184 134  84  70  70  27  12 

# 查看metadata
> head(liver[[]])
                    orig.ident nCount_RNA nFeature_RNA RNA_snn_res.1  seurat_clusters
AAACCTGAGCTACCTA SeuratProject       8116         1854             3                3
AAACCTGCAAGACACG SeuratProject       1748          913             4                4
AAACCTGCATGCCTTC SeuratProject       8838         2234             3                3
AAACCTGGTATCTGCA SeuratProject       1648          690             1                1
AAACGGGAGATATGGT SeuratProject       1583          499             0                0
AAACGGGTCCGCATAA SeuratProject       1492          705             1                1

> liver.metadata <- liver@meta.data
> liver.metadata$cell <- rownames(liver.metadata)
> head(liver.metadata)
                    orig.ident nCount_RNA nFeature_RNA RNA_snn_res.1 seurat_clusters             cell
AAACCTGAGCTACCTA SeuratProject       8116         1854             3               3 AAACCTGAGCTACCTA
AAACCTGCAAGACACG SeuratProject       1748          913             4               4 AAACCTGCAAGACACG
AAACCTGCATGCCTTC SeuratProject       8838         2234             3               3 AAACCTGCATGCCTTC
AAACCTGGTATCTGCA SeuratProject       1648          690             1               1 AAACCTGGTATCTGCA
AAACGGGAGATATGGT SeuratProject       1583          499             0               0 AAACGGGAGATATGGT
AAACGGGTCCGCATAA SeuratProject       1492          705             1               1 AAACGGGTCCGCATAA

> fin<-liver.metadata[,c(6,4)]
> head(fin)
                             cell RNA_snn_res.1
AAACCTGAGCTACCTA AAACCTGAGCTACCTA             3
AAACCTGCAAGACACG AAACCTGCAAGACACG             4
AAACCTGCATGCCTTC AAACCTGCATGCCTTC             3
AAACCTGGTATCTGCA AAACCTGGTATCTGCA             1
AAACGGGAGATATGGT AAACGGGAGATATGGT             0
AAACGGGTCCGCATAA AAACGGGTCCGCATAA             1

> write.csv(fin,file="Barcode_vs_ClusterID_of_liver_res1.csv", quote=FALSE, row.names=FALSE)

# 通过sftp 将"Barcode_vs_ClusterID_of_liver_res1.csv"文件下载到本地，用excel进行整合

========================================================================================================================================
## Lung
# Load the lung dataset
lung <- Read10X(data.dir = '/dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/droplet/Lung-10X_P7_9')
lung <- CreateSeuratObject(counts = lung)
lung
# An object of class Seurat 
# 23433 features across 1525 samples within 1 assay 
# Active assay: RNA (23433 features)

lung[["percent.mt"]] <- PercentageFeatureSet(object = lung, pattern = "^mt-")
pdf("lung_QC.pdf");VlnPlot(lung, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3);dev.off()
# Normalize the data
lung <- NormalizeData(lung)

# Identify highly variable features
lung <- FindVariableFeatures(lung, selection.method = "vst", nfeatures = 2000)

# Scale the data
all.genes <- rownames(lung)
lung <- ScaleData(lung, features = all.genes)

# Perform PCA on the scaled data
#lung <- RunPCA(lung, features = FindVariableFeatures(object = lung))
lung <- RunPCA(lung, npcs = 30, verbose = FALSE)

# Cluster the cells
lung <- RunUMAP(lung, reduction = "pca", dims = 1:20)
lung <- FindNeighbors(lung, reduction = "pca", dims = 1:20)
lung <- FindClusters(lung, resolution = 0.1); table(lung@active.ident)
#  0   1   2   3   4   5   6   7   8   9  10 
# 279 237 236 174 142 138 114  96  45  44  20  

# 查看metadata
head(lung[[]])
lung.metadata <- lung@meta.data
lung.metadata$cell <- rownames(lung.metadata)
head(lung.metadata)
fin<-lung.metadata[,c(10,9)]
head(fin)
                             cell RNA_snn_res.0.1
AAACCTGAGCCAGTAG AAACCTGAGCCAGTAG               1
AAACCTGGTAAACACA AAACCTGGTAAACACA               3
AAACCTGGTGAGTGAC AAACCTGGTGAGTGAC               2
AAACCTGGTTACGCGC AAACCTGGTTACGCGC               5
AAACCTGTCAGTACGT AAACCTGTCAGTACGT               7
AAACGGGCAATCTGCA AAACGGGCAATCTGCA               1

> write.csv(fin, file="Barcode_vs_ClusterID_of_lung_res0.1.csv", quote=FALSE, row.names=FALSE)


========================================================================================================================================
## Heart
# Load the heart dataset
heart <- Read10X(data.dir = '/dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/droplet/Heart-10X_P7_4')
heart <- CreateSeuratObject(counts = heart)
heart
# An object of class Seurat 
# 23433 features across 654 samples within 1 assay 
# Active assay: RNA (23433 features)

heart[["percent.mt"]] <- PercentageFeatureSet(object = heart, pattern = "^mt-")
pdf("heart_QC.pdf");VlnPlot(heart, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3);dev.off()
# Normalize the data
heart <- NormalizeData(heart)

# Identify highly variable features
heart <- FindVariableFeatures(heart, selection.method = "vst", nfeatures = 2000)

# Scale the data
all.genes <- rownames(heart)
heart <- ScaleData(heart, features = all.genes)

# Perform PCA on the scaled data
#heart <- RunPCA(heart, features = FindVariableFeatures(object = heart))
heart <- RunPCA(heart, npcs = 30, verbose = FALSE)

# Cluster the cells
heart <- RunUMAP(heart, reduction = "pca", dims = 1:20)
heart <- FindNeighbors(heart, reduction = "pca", dims = 1:20)
heart <- FindClusters(heart, resolution = 1.0); table(heart@active.ident)
#   0   1   2   3   4   5   6   7   8   9  10 
# 124  88  79  62  60  57  54  51  30  25  24 

# 查看metadata
head(heart[[]])
heart.metadata <- heart@meta.data
heart.metadata$cell <- rownames(heart.metadata)
head(heart.metadata)
fin<-heart.metadata[,c(7,5)]
head(fin)
                             cell RNA_snn_res.1
AAACCTGCACGACGAA AAACCTGCACGACGAA             0
AAACCTGGTTCCACGG AAACCTGGTTCCACGG             2
AAACCTGTCTCGATGA AAACCTGTCTCGATGA             0
AAACGGGAGCGCTCCA AAACGGGAGCGCTCCA             9
AAAGCAAAGTGAATTG AAAGCAAAGTGAATTG             7
AAAGCAACACCGAAAG AAAGCAACACCGAAAG             3

write.csv(fin,file="Barcode_vs_ClusterID_of_heart_res0.1.csv", quote=FALSE, row.names=FALSE)


========================================================================================================================================
## Kidney
# Load the kidney dataset
kidney <- Read10X(data.dir = '/dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/droplet/Kidney-10X_P7_5')
kidney <- CreateSeuratObject(counts = kidney)
kidney
# An object of class Seurat 
# 23433 features across 1264 samples within 1 assay thin 1 assay 
# Active assay: RNA (23433 features)

kidney[["percent.mt"]] <- PercentageFeatureSet(object = kidney, pattern = "^mt-")
pdf("kidney_QC.pdf");VlnPlot(kidney, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3);dev.off()
# Normalize the data
kidney <- NormalizeData(kidney)

# Identify highly variable features
kidney <- FindVariableFeatures(kidney, selection.method = "vst", nfeatures = 2000)

# Scale the data
all.genes <- rownames(kidney)
kidney <- ScaleData(kidney, features = all.genes)

# Perform PCA on the scaled data
#kidney <- RunPCA(kidney, features = FindVariableFeatures(object = kidney))
kidney <- RunPCA(kidney, npcs = 30, verbose = FALSE)

# Cluster the cells
kidney <- RunUMAP(kidney, reduction = "pca", dims = 1:20)
kidney <- FindNeighbors(kidney, reduction = "pca", dims = 1:20)
kidney <- FindClusters(kidney, resolution = 0.5); table(kidney@active.ident)
#  0   1   2   3   4   5   6   7   8   9  10 
# 347 278 225 110  88  85  43  27  24  22  15

# 查看metadata
head(kidney[[]])
kidney.metadata <- kidney@meta.data
kidney.metadata$cell <- rownames(kidney.metadata)
head(kidney.metadata)
                    orig.ident nCount_RNA nFeature_RNA RNA_snn_res.1 seurat_clusters             cell
AAACCTGAGCTACCTA SeuratProject       8116         1854             3               3 AAACCTGAGCTACCTA
AAACCTGCAAGACACG SeuratProject       1748          913             4               4 AAACCTGCAAGACACG
AAACCTGCATGCCTTC SeuratProject       8838         2234             3               3 AAACCTGCATGCCTTC
AAACCTGGTATCTGCA SeuratProject       1648          690             1               1 AAACCTGGTATCTGCA
AAACGGGAGATATGGT SeuratProject       1583          499             0               0 AAACGGGAGATATGGT
AAACGGGTCCGCATAA SeuratProject       1492          705             1               1 AAACGGGTCCGCATAA

fin<-kidney.metadata[,c(8,7)]
head(fin)
                             cell RNA_snn_res.0.5
AAACCTGCATGTTGAC AAACCTGCATGTTGAC               0
AAACGGGGTGCGGTAA AAACGGGGTGCGGTAA               2
AAACGGGTCATGCATG AAACGGGTCATGCATG               9
AAAGATGAGGACCACA AAAGATGAGGACCACA               1
AAAGATGTCAACACAC AAAGATGTCAACACAC               4
AAAGCAAAGCCACGTC AAAGCAAAGCCACGTC               1

write.csv(fin,file="Barcode_vs_ClusterID_of_kidney_res0.5.csv", quote=FALSE, row.names=FALSE)

# 通过sftp 将"Barcode_vs_ClusterID_of_heart_res0.1.csv"文件下载到本地，用excel可以打开
get /dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/seurat_liver/Barcode_vs_ClusterID_of_kidney_res0.5.csv ./Desktop
get /dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/seurat_liver/Barcode_vs_ClusterID_of_heart_res0.1.csv ./Desktop
get /dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/seurat_liver/Barcode_vs_ClusterID_of_lung_res0.1.csv ./Desktop
get /dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/seurat_liver/Barcode_vs_ClusterID_of_liver_res1.csv ./Desktop



# 通过sftp 将"Heart-10X_P7_4.csv"文件下载到本地，用excel可以打开
get /dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/artifical_matrix/Heart-10X_P7_4.csv ./Desktop
get /dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/artifical_matrix/Liver/Liver-10X_P4_2.csv ./Desktop 
get /dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/artifical_matrix/Lung-10X_P7_9.csv ./Desktop
get /dellfsqd2/ST_LBI/USER/liangxue/project/10x_transcriptomics/Mouse/07.download_data/20200113/artifical_matrix/Kidney-10X_P7_5.csv ./Desktop






