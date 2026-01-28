# 设置详细输出
options(echo = TRUE)
# 记录开始时间
#cat("Script started at:", Sys.time(), "\n")

#设置镜像
local({r<-getOption("repos")  
r["CRAN"]<-"https://mirrors.tuna.tsinghua.edu.cn/CRAN/"   
options(repos = r)})

#检查是否安装meta包，若未安装，则自动安装
p="meta"
if(!suppressWarnings(suppressMessages(require(p,character.only = TRUE,quietly = TRUE,warn.conflicts = FALSE)))){
  install.packages(p,warn.conflicts = FALSE)
  suppressWarnings(suppressMessages(library(p,character.only = TRUE,quietly = TRUE,warn.conflicts = FALSE)))
}
#检查是否安装optparse包，若未安装，则自动安装
p="optparse"
if(!suppressWarnings(suppressMessages(require(p,character.only = TRUE,quietly = TRUE,warn.conflicts = FALSE)))){
  install.packages(p,warn.conflicts = FALSE)
  suppressWarnings(suppressMessages(library(p,character.only = TRUE,quietly = TRUE,warn.conflicts = FALSE)))
}

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
 	  #-m/--method：树状图绘制方法，默认倒方差法（I）
 	  make_option(c("-m", "--method"), type = "character", default = "I", 
				  help = "树状图绘制方法 [默认: %default]", metavar = "METHOD"),
		#-p/--pool：敏感性分析计算模型，默认随机效应模型（random）
 	  make_option(c("-p", "--pool"), type = "character", default = "random", 
				  help = "敏感性分析计算模型 [默认: %default]", metavar = "METHOD"), 	  
 	  #-o/--outdir：输出目录，默认当前工作目录
 	  make_option(c("-o", "--outdir"), type = "character", default = getwd(),
				  help = "输出目录 [默认: 当前目录]", metavar = "DIR")
)
#Rscript galbraith_text.R -l mian.csv -g subgroup_children_adults.csv -s severity.csv

# 解析参数
opt <- parse_args(OptionParser(option_list = option_list))

# 检查必需参数
if (is.null(opt$maindata) || is.null(opt$subdata) || is.null(opt$severitydata)) {
  stop("必须提供 -l/--maindata 和 -g/--subdata 参数 以及 -s/--severitydata 参数")
}

# 创建输出目录（如果不存在）
if (!dir.exists(opt$outdir)) {
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
}


#////////////////////////////
#数据可视化绘制

# 包的加载
package_list <- "meta"

Main<-read.csv(opt$maindata,header=TRUE,dec=".")
Sub<-read.csv(opt$subdata,header=TRUE,dec=".")
Severity<-read.csv(opt$severitydata,header=TRUE,dec=".")

# Dominant model galbraith 显性模型甘氏图(定位离群值)

event.e<-Main$AA_AG_case
n.e<-Main$AA_AG_case+Main$GG_case
event.c<-Main$AA_AG_control
n.c<-Main$AA_AG_control+Main$GG_control

meta1<-metabin(event.e,n.e,event.c,n.c,data=Main,
	method=opt$method,sm="OR", #method = "MH"：指定使用 Mantel-Haenszel法 计算合并效应量
										#特点：适用于研究间效应量同质性较好的情况
										#优势：在小样本或零单元格（zero-cell）情况下表现稳健
										#其他可选方法： 1."Inverse"：逆方差法（默认） 2."Peto"：Peto法，适用于罕见事件（OR接近1） 3."GLMM"：广义线性混合模型
	common=TRUE,random=TRUE,
	studlab=paste(Main$Author,Main$Year))

FileName<-paste("Dominant model galbraith_text.PDF")
pdf(file=paste(opt$outdir,"/","Dominant_Model","/",FileName,sep = ""),height=8,width=8)

#radial(meta1,level=0.95) #设置置信区间为95%
radial(meta1,text=Main$Author,level=0.95)

dev.off()

event.e<-Severity$AA_AG_severe
n.e<-Severity$AA_AG_severe+Severity$GG_severe
event.c<-Severity$AA_AG_not_severe
n.c<-Severity$GG_not_severe+Severity$AA_AG_not_severe

meta1<-metabin(event.e,n.e,event.c,n.c,data=Severity,
	method=opt$method,sm="OR", #method = "MH"：指定使用 Mantel-Haenszel法 计算合并效应量
										#特点：适用于研究间效应量同质性较好的情况
										#优势：在小样本或零单元格（zero-cell）情况下表现稳健
										#其他可选方法： 1."Inverse"：逆方差法（默认） 2."Peto"：Peto法，适用于罕见事件（OR接近1） 3."GLMM"：广义线性混合模型
	common=TRUE,random=TRUE,
	studlab=paste(Severity$Author,Severity$Year))

FileName<-paste("Dominant model galbraith_severity_text.PDF")
pdf(file=paste(opt$outdir,"/","Dominant_Model","/",FileName,sep = ""),height=8,width=8)

#radial(meta1,level=0.95) #设置置信区间为95%
radial(meta1,text=Severity$Author,level=0.95)

dev.off()

# Heterozygote model galbraith 杂合子模型甘氏图(定位离群值)

event.e<-Main$AG_case
n.e<-Main$AG_case+Main$AA_case+Main$GG_case
event.c<-Main$AG_control
n.c<-Main$AG_control+Main$AA_control+Main$GG_control

meta1<-metabin(event.e,n.e,event.c,n.c,data=Main,
	method=opt$method,sm="OR", #method = "MH"：指定使用 Mantel-Haenszel法 计算合并效应量
										#特点：适用于研究间效应量同质性较好的情况
										#优势：在小样本或零单元格（zero-cell）情况下表现稳健
										#其他可选方法： 1."Inverse"：逆方差法（默认） 2."Peto"：Peto法，适用于罕见事件（OR接近1） 3."GLMM"：广义线性混合模型
	common=TRUE,random=TRUE,
	studlab=paste(Main$Author,Main$Year))

FileName<-paste("Heterozygote model galbraith_text.PDF")
pdf(file=paste(opt$outdir,"/","Heterozygote_Model","/",FileName,sep = ""),height=8,width=8)

#radial(meta1,level=0.95) #设置置信区间为95%
radial(meta1,text=Main$Author,level=0.95)

dev.off()

event.e<-Severity$AG_severe
n.e<-Severity$AG_severe+Severity$AA_severe+Severity$GG_severe
event.c<-Severity$AG_not_severe
n.c<-Severity$AA_not_severe+Severity$AG_not_severe+Severity$GG_not_severe

meta1<-metabin(event.e,n.e,event.c,n.c,data=Severity,
	method=opt$method,sm="OR", #method = "MH"：指定使用 Mantel-Haenszel法 计算合并效应量
										#特点：适用于研究间效应量同质性较好的情况
										#优势：在小样本或零单元格（zero-cell）情况下表现稳健
										#其他可选方法： 1."Inverse"：逆方差法（默认） 2."Peto"：Peto法，适用于罕见事件（OR接近1） 3."GLMM"：广义线性混合模型
	common=TRUE,random=TRUE,
	studlab=paste(Severity$Author,Severity$Year))

FileName<-paste("Heterozygote model galbraith_severity_text.PDF")
pdf(file=paste(opt$outdir,"/","Heterozygote_Model","/",FileName,sep = ""),height=8,width=8)

#radial(meta1,level=0.95) #设置置信区间为95%
radial(meta1,text=Severity$Author,level=0.95)

dev.off()

# Homozygote model galbraith 纯合子模型甘氏图(定位离群值)

event.e<-Main$AA_case
n.e<-Main$AA_case+Main$GG_case
event.c<-Main$AA_control
n.c<-Main$AA_control+Main$GG_control

meta1<-metabin(event.e,n.e,event.c,n.c,data=Main,
	method=opt$method,sm="OR", #method = "MH"：指定使用 Mantel-Haenszel法 计算合并效应量
										#特点：适用于研究间效应量同质性较好的情况
										#优势：在小样本或零单元格（zero-cell）情况下表现稳健
										#其他可选方法： 1."Inverse"：逆方差法（默认） 2."Peto"：Peto法，适用于罕见事件（OR接近1） 3."GLMM"：广义线性混合模型
	common=TRUE,random=TRUE,
	studlab=paste(Main$Author,Main$Year))

FileName<-paste("Homozygote model galbraith_text.PDF")
pdf(file=paste(opt$outdir,"/","Homozygote_Model","/",FileName,sep = ""),height=8,width=8)

#radial(meta1,level=0.95) #设置置信区间为95%
radial(meta1,text=Main$Author,level=0.95)

dev.off()

event.e<-Severity$AA_severe
n.e<-Severity$AA_severe+Severity$GG_severe
event.c<-Severity$AA_not_severe
n.c<-Severity$AA_not_severe+Severity$GG_not_severe

meta1<-metabin(event.e,n.e,event.c,n.c,data=Severity,
	method=opt$method,sm="OR", #method = "MH"：指定使用 Mantel-Haenszel法 计算合并效应量
										#特点：适用于研究间效应量同质性较好的情况
										#优势：在小样本或零单元格（zero-cell）情况下表现稳健
										#其他可选方法： 1."Inverse"：逆方差法（默认） 2."Peto"：Peto法，适用于罕见事件（OR接近1） 3."GLMM"：广义线性混合模型
	common=TRUE,random=TRUE,
	studlab=paste(Severity$Author,Severity$Year))

FileName<-paste("Homozygote model galbraith_severity_text.PDF")
pdf(file=paste(opt$outdir,"/","Homozygote_Model","/",FileName,sep = ""),height=8,width=8)

#radial(meta1,level=0.95) #设置置信区间为95%
radial(meta1,text=Severity$Author,level=0.95)

dev.off()

# Recessive model galbraith 隐性模型甘氏图(定位离群值)

event.e<-Main$AA_case
n.e<-Main$AA_case+Main$AG_GG_case
event.c<-Main$AA_control
n.c<-Main$AA_control+Main$AG_GG_control

meta1<-metabin(event.e,n.e,event.c,n.c,data=Main,
	method=opt$method,sm="OR", #method = "MH"：指定使用 Mantel-Haenszel法 计算合并效应量
										#特点：适用于研究间效应量同质性较好的情况
										#优势：在小样本或零单元格（zero-cell）情况下表现稳健
										#其他可选方法： 1."Inverse"：逆方差法（默认） 2."Peto"：Peto法，适用于罕见事件（OR接近1） 3."GLMM"：广义线性混合模型
	common=TRUE,random=TRUE,
	studlab=paste(Main$Author,Main$Year))

FileName<-paste("Recessive model galbraith_text.PDF")
pdf(file=paste(opt$outdir,"/","Recessive_Model","/",FileName,sep = ""),height=8,width=8)

#radial(meta1,level=0.95) #设置置信区间为95%
radial(meta1,text=Main$Author,level=0.95)

dev.off()

event.e<-Severity$AA_severe
n.e<-Severity$AA_severe+Severity$AG_GG_severe
event.c<-Severity$AA_not_severe
n.c<-Severity$AA_not_severe+Severity$AG_GG_not_severe

meta1<-metabin(event.e,n.e,event.c,n.c,data=Severity,
	method=opt$method,sm="OR", #method = "MH"：指定使用 Mantel-Haenszel法 计算合并效应量
										#特点：适用于研究间效应量同质性较好的情况
										#优势：在小样本或零单元格（zero-cell）情况下表现稳健
										#其他可选方法： 1."Inverse"：逆方差法（默认） 2."Peto"：Peto法，适用于罕见事件（OR接近1） 3."GLMM"：广义线性混合模型
	common=TRUE,random=TRUE,
	studlab=paste(Severity$Author,Severity$Year))

FileName<-paste("Recessive model galbraith_severity_text.PDF")
pdf(file=paste(opt$outdir,"/","Recessive_Model","/",FileName,sep = ""),height=8,width=8)

#radial(meta1,level=0.95) #设置置信区间为95%
radial(meta1,text=Severity$Author,level=0.95)

dev.off()

# Allele model galbraith 等位基因模型甘氏图(定位离群值)

event.e<-Main$A_case
n.e<-Main$A_case+Main$G_case
event.c<-Main$A_control
n.c<-Main$A_control+Main$G_control

meta1<-metabin(event.e,n.e,event.c,n.c,data=Main,
	method=opt$method,sm="OR", #method = "MH"：指定使用 Mantel-Haenszel法 计算合并效应量
										#特点：适用于研究间效应量同质性较好的情况
										#优势：在小样本或零单元格（zero-cell）情况下表现稳健
										#其他可选方法： 1."Inverse"：逆方差法（默认） 2."Peto"：Peto法，适用于罕见事件（OR接近1） 3."GLMM"：广义线性混合模型
	common=TRUE,random=TRUE,
	studlab=paste(Main$Author,Main$Year))

FileName<-paste("Allele model galbraith_text.PDF")
pdf(file=paste(opt$outdir,"/","Allele_Model","/",FileName,sep = ""),height=8,width=8)

#radial(meta1,level=0.95) #设置置信区间为95%
radial(meta1,text=Main$Author,level=0.95)

dev.off()

event.e<-Severity$A_severe
n.e<-Severity$A_severe+Severity$G_severe
event.c<-Severity$A_not_severe
n.c<-Severity$A_not_severe+Severity$G_not_severe

meta1<-metabin(event.e,n.e,event.c,n.c,data=Severity,
	method=opt$method,sm="OR", #method = "MH"：指定使用 Mantel-Haenszel法 计算合并效应量
										#特点：适用于研究间效应量同质性较好的情况
										#优势：在小样本或零单元格（zero-cell）情况下表现稳健
										#其他可选方法： 1."Inverse"：逆方差法（默认） 2."Peto"：Peto法，适用于罕见事件（OR接近1） 3."GLMM"：广义线性混合模型
	common=TRUE,random=TRUE,
	studlab=paste(Severity$Author,Severity$Year))

FileName<-paste("Allele model galbraith_severity_text.PDF")
pdf(file=paste(opt$outdir,"/","Allele_Model","/",FileName,sep = ""),height=8,width=8)

#radial(meta1,level=0.95) #设置置信区间为95%
radial(meta1,text=Severity$Author,level=0.95)

dev.off()
