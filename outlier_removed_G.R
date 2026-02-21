# 设置详细输出
#options(echo = TRUE)

# 使用 optparse 解码命令行参数
if (!require("optparse")) {
  install.packages("optparse", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  library(optparse)
}

# 定义命令行参数
	option_list <- list(
 	  #-l/--data：主要数据文件路径（必需）
 	  make_option(c("-l", "--maindata"), type = "character", 
				  help = "主要数据文件[必需]", metavar = "FILE"),
 	  #-g/--subdata：亚组数据文件路径（必需）
 	  make_option(c("-g", "--subdata"), type = "character",
				  help = "亚组数据文件[必需]", metavar = "FILE"),
 	  #-s/--severitydata：严重程度数据文件路径（必需）
 	  make_option(c("-s", "--severitydata"), type = "character",
				  help = "严重程度数据文件[必需]", metavar = "FILE"),
 	  make_option(c("-r", "--rdata"), type = "character",
				  help = "严重程度数据文件[必需]", metavar = "FILE"),
 	  make_option(c("-a", "--adata"), type = "character",
				  help = "严重程度数据文件[必需]", metavar = "FILE"),
 	  #-m/--method：树状图绘制方法，默认MH法
 	  make_option(c("-m", "--method"), type = "character", default = "MH", 
				  help = "树状图绘制方法 [默认: %default]", metavar = "METHOD"),
		#-p/--pool：敏感性分析计算模型，默认随机效应模型（random）
 	  make_option(c("-p", "--pool"), type = "character", default = "random", 
				  help = "敏感性分析计算模型 [默认: %default]", metavar = "METHOD"), 	  
 	  #-o/--outdir：输出目录，默认当前工作目录
 	  make_option(c("-o", "--outdir"), type = "character", default = getwd(),
				  help = "输出目录 [默认: 当前目录]", metavar = "DIR")
)
#Rscript r.meta1_reverse.R -l main_D.csv -g main_He.csv -s main_Ho.csv -r main_R.csv -a main_A.csv

# 解析参数
opt <- parse_args(OptionParser(option_list = option_list))

# 创建输出目录（如果不存在）
if (!dir.exists(opt$outdir)) {
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
}

library("meta")
## Dominant model forest 显性模型
cat("--- Dominant model forest 显性模型 ---\n")
Main<-read.csv(opt$maindata,header=TRUE,dec=".")

event.e<-Main$AG_GG_case
n.e<-Main$AA_AG_case+Main$GG_case
event.c<-Main$AG_GG_control
n.c<-Main$AA_AG_control+Main$GG_control

meta1<-metabin(event.e,n.e,event.c,n.c,data=Main,
	method=opt$method,sm="OR", #method = "MH"：指定使用 Mantel-Haenszel法 计算合并效应量
										#特点：适用于研究间效应量同质性较好的情况
										#优势：在小样本或零单元格（zero-cell）情况下表现稳健
										#其他可选方法： 1."Inverse"：逆方差法（默认） 2."Peto"：Peto法，适用于罕见事件（OR接近1） 3."GLMM"：广义线性混合模型
	common=TRUE,random=TRUE,
	studlab=paste(Main$Author,Main$Year))

meta1

## Heterozygote model forest 杂合子模型
cat("--- Heterozygote model forest 杂合子模型 ---\n")
Main<-read.csv(opt$subdata,header=TRUE,dec=".")

event.e<-Main$AG_case
n.e<-Main$AG_case+Main$AA_case
event.c<-Main$AG_control
n.c<-Main$AG_control+Main$AA_control

meta1<-metabin(event.e,n.e,event.c,n.c,data=Main,
	method=opt$method,sm="OR", #method = "MH"：指定使用 Mantel-Haenszel法 计算合并效应量
										#特点：适用于研究间效应量同质性较好的情况
										#优势：在小样本或零单元格（zero-cell）情况下表现稳健
										#其他可选方法： 1."Inverse"：逆方差法（默认） 2."Peto"：Peto法，适用于罕见事件（OR接近1） 3."GLMM"：广义线性混合模型
	common=TRUE,random=TRUE,
	studlab=paste(Main$Author,Main$Year))

meta1

## Homozygote model forest 纯合子模型
cat("--- Homozygote model forest 纯合子模型 ---\n")
Main<-read.csv(opt$severitydata,header=TRUE,dec=".")

event.e<-Main$GG_case
n.e<-Main$AA_case+Main$GG_case
event.c<-Main$GG_control
n.c<-Main$AA_control+Main$GG_control

meta1<-metabin(event.e,n.e,event.c,n.c,data=Main,
	method=opt$method,sm="OR", #method = "MH"：指定使用 Mantel-Haenszel法 计算合并效应量
										#特点：适用于研究间效应量同质性较好的情况
										#优势：在小样本或零单元格（zero-cell）情况下表现稳健
										#其他可选方法： 1."Inverse"：逆方差法（默认） 2."Peto"：Peto法，适用于罕见事件（OR接近1） 3."GLMM"：广义线性混合模型
	common=TRUE,random=TRUE,
	studlab=paste(Main$Author,Main$Year))

meta1

## Recessive model forest 隐性模型
cat("--- Recessive model forest 隐性模型 ---\n")
Main<-read.csv(opt$rdata,header=TRUE,dec=".")

event.e<-Main$GG_case
n.e<-Main$AA_case+Main$AG_GG_case
event.c<-Main$GG_control
n.c<-Main$AA_control+Main$AG_GG_control

meta1<-metabin(event.e,n.e,event.c,n.c,data=Main,
	method=opt$method,sm="OR", #method = "MH"：指定使用 Mantel-Haenszel法 计算合并效应量
										#特点：适用于研究间效应量同质性较好的情况
										#优势：在小样本或零单元格（zero-cell）情况下表现稳健
										#其他可选方法： 1."Inverse"：逆方差法（默认） 2."Peto"：Peto法，适用于罕见事件（OR接近1） 3."GLMM"：广义线性混合模型
	common=TRUE,random=TRUE,
	studlab=paste(Main$Author,Main$Year))

meta1

##Allele model 等位基因模型
cat("--- Allele model 等位基因模型 ---\n")
Main<-read.csv(opt$adata,header=TRUE,dec=".")

event.e<-Main$G_case
n.e<-Main$A_case+Main$G_case
event.c<-Main$G_control
n.c<-Main$A_control+Main$G_control

meta1<-metabin(event.e,n.e,event.c,n.c,data=Main,
	method=opt$method,sm="OR", #method = "MH"：指定使用 Mantel-Haenszel法 计算合并效应量
										#特点：适用于研究间效应量同质性较好的情况
										#优势：在小样本或零单元格（zero-cell）情况下表现稳健
										#其他可选方法： 1."Inverse"：逆方差法（默认） 2."Peto"：Peto法，适用于罕见事件（OR接近1） 3."GLMM"：广义线性混合模型
	common=TRUE,random=TRUE,
	studlab=paste(Main$Author,Main$Year))

meta1