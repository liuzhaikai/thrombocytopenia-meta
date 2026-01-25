# @Author: Zhaikai Liu
# @Email:  lzk_kklk@outlook.com
# @Timestamp for Creation: 2026-01-24 22:20:00
# @Last Modified by: Zhaikai Liu
# @Last Modified time: 2026-01-25 10：01：30

# IL-17F rs763780  显性：T(A) 隐性：C(G)

#注：暴露组标签设置为"Exposure"；对照组标签设置为"Control"

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
	#-l/--data：主要数据文件路径（必需）
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
#Rscript thrombocytopenia_meta_1.R -l mian.csv -g subgroup_children_adults.csv -s severity.csv
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

## Dominant model forest 显性模型

Main<-read.csv(opt$maindata,header=TRUE,dec=".")
Sub<-read.csv(opt$subdata,header=TRUE,dec=".")
Severity<-read.csv(opt$severitydata,header=TRUE,dec=".")

event.e<-Main$AG_GG_case
n.e<-Main$AG_GG_case+Main$AG_case
event.c<-Main$AG_GG_control
n.c<-Main$AG_GG_control+Main$AA_control

meta1<-metabin(event.e,n.e,event.c,n.c,data=Main,
	method=opt$method,sm="OR", #method = "MH"：指定使用 Mantel-Haenszel法 计算合并效应量
										#特点：适用于研究间效应量同质性较好的情况
										#优势：在小样本或零单元格（zero-cell）情况下表现稳健
										#其他可选方法： 1."Inverse"：逆方差法（默认） 2."Peto"：Peto法，适用于罕见事件（OR接近1） 3."GLMM"：广义线性混合模型
	common=TRUE,random=TRUE,
	studlab=paste(Main$Author,Main$Year))

FileName<-paste("Dominant model forest.PDF")
pdf(file=paste(opt$outdir,"/","Dominant_Model","/",FileName,sep = ""),height=4.5,width=12)

forest(meta1,lab.e="Exposure",lab.c="Control")#lab.e:暴露组设置自定义标签（Exposure） lab.c:为对照组设置自定义标签（Control）

dev.off()

# Dominant model sensetivity_analysis

#敏感性分析，注意是用固定效应模型还是随机效应模型
FileName<-paste("Dominant model sensetivity_analysis_",opt$pool,".PDF")
pdf(file=paste(opt$outdir,"/","Dominant_Model","/",FileName,sep = ""),height=3.5,width=10)

forest(metainf(meta1,pooled=opt$pool),lab.e="Exposure",lab.c="Control")

dev.off()

# Dominant model funnel 显性模型漏斗图

FileName<-paste("Dominant model funnel.PDF")
pdf(file=paste(opt$outdir,"/","Dominant_Model","/",FileName,sep = ""),height=8,width=8)

funnel(meta1)

dev.off()

# Dominant model galbraith 显性模型甘氏图(定位离群值)

FileName<-paste("Dominant model galbraith.PDF")
pdf(file=paste(opt$outdir,"/","Dominant_Model","/",FileName,sep = ""),height=8,width=8)

# radial(meta1,level=0.95) #设置置信区间为95%
radial(meta1,text=Main$Author,level=0.95)

dev.off()

# Dominant model subgroup_analysis_Country forest

FileName<-paste("Dominant model subgroup_analysis_Country.PDF")
pdf(file=paste(opt$outdir,"/","Dominant_Model","/",FileName,sep = ""),height=9.5,width=14)

meta1<-metabin(event.e,n.e,event.c,n.c,data=Main,
	method=opt$method,sm="OR",byvar=Main$Country,
	common=TRUE,random=TRUE,
	studlab=paste(Main$Author,Main$Year))

forest(meta1,lab.e="Exposure",lab.c="Control")

dev.off()

# Dominant model subgroup_analysis_Sample_age forest

event.e<-Sub$AG_GG_case
n.e<-Sub$AG_GG_case+Sub$AG_case
event.c<-Sub$AG_GG_control
n.c<-Sub$AG_GG_control+Sub$AA_control

FileName<-paste("Dominant model subgroup_analysis_Sample_age.PDF")
pdf(file=paste(opt$outdir,"/","Dominant_Model","/",FileName,sep = ""),height=9.5,width=14)

meta1<-metabin(event.e,n.e,event.c,n.c,data=Sub,
	method=opt$method,sm="OR",byvar=Sub$Sample_age,
	common=TRUE,random=TRUE,
	studlab=paste(Sub$Author,Sub$Year))

forest(meta1,lab.e="Exposure",lab.c="Control")

dev.off()

# Dominant model forest_severity

event.e<-Severity$AG_GG_not_severe
n.e<-Severity$AG_GG_not_severe+Severity$AG_not_severe
event.c<-Severity$AG_GG_severe
n.c<-Severity$AG_GG_severe+Severity$AA_severe

meta1<-metabin(event.e,n.e,event.c,n.c,data=Severity,
	method=opt$method,sm="OR", #method = "MH"：指定使用 Mantel-Haenszel法 计算合并效应量
										#特点：适用于研究间效应量同质性较好的情况
										#优势：在小样本或零单元格（zero-cell）情况下表现稳健
										#其他可选方法： 1."Inverse"：逆方差法（默认） 2."Peto"：Peto法，适用于罕见事件（OR接近1） 3."GLMM"：广义线性混合模型
	common=TRUE,random=TRUE,
	studlab=paste(Severity$Author,Severity$Year))

FileName<-paste("Dominant model forest_severity.PDF")
pdf(file=paste(opt$outdir,"/","Dominant_Model","/",FileName,sep = ""),height=4.5,width=12)

forest(meta1,lab.e="Not_Severe",lab.c="Severe")

dev.off()

# Dominant model sensetivity_analysis_severity

FileName<-paste("Dominant model sensetivity_analysis_",opt$pool,"_severity.PDF")
pdf(file=paste(opt$outdir,"/","Dominant_Model","/",FileName,sep = ""),height=3.5,width=10)

forest(metainf(meta1,pooled=opt$pool),lab.e="Not_Severe",lab.c="Severe")

dev.off()

# Dominant model funnel_severity 

FileName<-paste("Dominant model funnel_severity.PDF")
pdf(file=paste(opt$outdir,"/","Dominant_Model","/",FileName,sep = ""),height=8,width=8)

funnel(meta1)

dev.off()

# Dominant model galbraith_severity

FileName<-paste("Dominant model galbraith_severity.PDF")
pdf(file=paste(opt$outdir,"/","Dominant_Model","/",FileName,sep = ""),height=8,width=8)

# radial(meta1,level=0.95) #设置置信区间为95%
radial(meta1,text=Severity$Author,level=0.95)

dev.off()


## Heterozygote model forest 杂合子模型

event.e<-Main$AG_case
n.e<-Main$AA_AG_case
event.c<-Main$AG_control
n.c<-Main$AA_AG_control

meta1<-metabin(event.e,n.e,event.c,n.c,data=Main,
	method=opt$method,sm="OR", #method = "MH"：指定使用 Mantel-Haenszel法 计算合并效应量
										#特点：适用于研究间效应量同质性较好的情况
										#优势：在小样本或零单元格（zero-cell）情况下表现稳健
										#其他可选方法： 1."Inverse"：逆方差法（默认） 2."Peto"：Peto法，适用于罕见事件（OR接近1） 3."GLMM"：广义线性混合模型
	common=TRUE,random=TRUE,
	studlab=paste(Main$Author,Main$Year))

FileName<-paste("Heterozygote model forest.PDF")
pdf(file=paste(opt$outdir,"/","Heterozygote_Model","/",FileName,sep = ""),height=4.5,width=12)

forest(meta1,lab.e="Exposure",lab.c="Control")#lab.e:暴露组设置自定义标签（Exposure） lab.c:为对照组设置自定义标签（Control）

dev.off()

# Heterozygote model sensetivity_analysis

#敏感性分析，注意是用固定效应模型还是随机效应模型
FileName<-paste("Heterozygote model sensetivity_analysis_",opt$pool,".PDF")
pdf(file=paste(opt$outdir,"/","Heterozygote_Model","/",FileName,sep = ""),height=3.5,width=10)

forest(metainf(meta1,pooled=opt$pool),lab.e="Exposure",lab.c="Control")

dev.off()

# Heterozygote model funnel 杂合子模型漏斗图

FileName<-paste("Heterozygote model funnel.PDF")
pdf(file=paste(opt$outdir,"/","Heterozygote_Model","/",FileName,sep = ""),height=8,width=8)

funnel(meta1)

dev.off()

# Heterozygote model galbraith 杂合子模型甘氏图(定位离群值)

FileName<-paste("Heterozygote model galbraith.PDF")
pdf(file=paste(opt$outdir,"/","Heterozygote_Model","/",FileName,sep = ""),height=8,width=8)

# radial(meta1,level=0.95) #设置置信区间为95%
radial(meta1,text=Main$Author,level=0.95)

dev.off()

# Heterozygote model subgroup_analysis_Country forest

FileName<-paste("Heterozygote model subgroup_analysis_Country.PDF")
pdf(file=paste(opt$outdir,"/","Heterozygote_Model","/",FileName,sep = ""),height=9.5,width=14)

meta1<-metabin(event.e,n.e,event.c,n.c,data=Main,
	method=opt$method,sm="OR",byvar=Main$Country,
	common=TRUE,random=TRUE,
	studlab=paste(Main$Author,Main$Year))

forest(meta1,lab.e="Exposure",lab.c="Control")

dev.off()

# Heterozygote model subgroup_analysis_Sample_age forest

event.e<-Sub$AG_case
n.e<-Sub$AA_AG_case
event.c<-Sub$AG_control
n.c<-Sub$AA_AG_control

FileName<-paste("Heterozygote model subgroup_analysis_Sample_age.PDF")
pdf(file=paste(opt$outdir,"/","Heterozygote_Model","/",FileName,sep = ""),height=9.5,width=14)

meta1<-metabin(event.e,n.e,event.c,n.c,data=Sub,
	method=opt$method,sm="OR",byvar=Sub$Sample_age,
	common=TRUE,random=TRUE,
	studlab=paste(Sub$Author,Sub$Year))

forest(meta1,lab.e="Exposure",lab.c="Control")

dev.off()

# Heterozygote model forest_severity

event.e<-Severity$AG_not_severe
n.e<-Severity$AA_AG_not_severe
event.c<-Severity$AG_severe
n.c<-Severity$AA_AG_severe

meta1<-metabin(event.e,n.e,event.c,n.c,data=Severity,
	method=opt$method,sm="OR", #method = "MH"：指定使用 Mantel-Haenszel法 计算合并效应量
										#特点：适用于研究间效应量同质性较好的情况
										#优势：在小样本或零单元格（zero-cell）情况下表现稳健
										#其他可选方法： 1."Inverse"：逆方差法（默认） 2."Peto"：Peto法，适用于罕见事件（OR接近1） 3."GLMM"：广义线性混合模型
	common=TRUE,random=TRUE,
	studlab=paste(Severity$Author,Severity$Year))

FileName<-paste("Heterozygote model forest_severity.PDF")
pdf(file=paste(opt$outdir,"/","Heterozygote_Model","/",FileName,sep = ""),height=4.5,width=12)

forest(meta1,lab.e="Not_Severe",lab.c="Severe")

dev.off()

# Heterozygote model sensetivity_analysis_severity

FileName<-paste("Heterozygote model sensetivity_analysis_",opt$pool,"_severity.PDF")
pdf(file=paste(opt$outdir,"/","Heterozygote_Model","/",FileName,sep = ""),height=3.5,width=10)

forest(metainf(meta1,pooled=opt$pool),lab.e="Not_Severe",lab.c="Severe")

dev.off()

# Heterozygote model funnel_severity 

FileName<-paste("Heterozygote model funnel_severity.PDF")
pdf(file=paste(opt$outdir,"/","Heterozygote_Model","/",FileName,sep = ""),height=8,width=8)

funnel(meta1)

dev.off()

# Heterozygote model galbraith_severity

FileName<-paste("Heterozygote model galbraith_severity.PDF")
pdf(file=paste(opt$outdir,"/","Heterozygote_Model","/",FileName,sep = ""),height=8,width=8)

# radial(meta1,level=0.95) #设置置信区间为95%
radial(meta1,text=Severity$Author,level=0.95)

dev.off()


## Homozygote model forest 纯合子模型

event.e<-Main$GG_case
n.e<-Main$GG_case+Main$AA_case
event.c<-Main$GG_control
n.c<-Main$GG_control+Main$AA_control

meta1<-metabin(event.e,n.e,event.c,n.c,data=Main,
	method=opt$method,sm="OR", #method = "MH"：指定使用 Mantel-Haenszel法 计算合并效应量
										#特点：适用于研究间效应量同质性较好的情况
										#优势：在小样本或零单元格（zero-cell）情况下表现稳健
										#其他可选方法： 1."Inverse"：逆方差法（默认） 2."Peto"：Peto法，适用于罕见事件（OR接近1） 3."GLMM"：广义线性混合模型
	common=TRUE,random=TRUE,
	studlab=paste(Main$Author,Main$Year))

FileName<-paste("Homozygote model forest.PDF")
pdf(file=paste(opt$outdir,"/","Homozygote_Model","/",FileName,sep = ""),height=4.5,width=12)

forest(meta1,lab.e="Exposure",lab.c="Control")#lab.e:暴露组设置自定义标签（Exposure） lab.c:为对照组设置自定义标签（Control）

dev.off()

# Homozygote model sensetivity_analysis

#敏感性分析，注意是用固定效应模型还是随机效应模型
FileName<-paste("Homozygote model sensetivity_analysis_",opt$pool,".PDF")
pdf(file=paste(opt$outdir,"/","Homozygote_Model","/",FileName,sep = ""),height=3.5,width=10)

forest(metainf(meta1,pooled=opt$pool),lab.e="Exposure",lab.c="Control")

dev.off()

# Homozygote model funnel 纯合子模型漏斗图

FileName<-paste("Homozygote model funnel.PDF")
pdf(file=paste(opt$outdir,"/","Homozygote_Model","/",FileName,sep = ""),height=8,width=8)

funnel(meta1)

dev.off()

# Homozygote model galbraith 纯合子模型甘氏图(定位离群值)

FileName<-paste("Homozygote model galbraith.PDF")
pdf(file=paste(opt$outdir,"/","Homozygote_Model","/",FileName,sep = ""),height=8,width=8)

# radial(meta1,level=0.95) #设置置信区间为95%
radial(meta1,text=Main$Author,level=0.95)

dev.off()

# Homozygote model subgroup_analysis_Country forest

FileName<-paste("Homozygote model subgroup_analysis_Country.PDF")
pdf(file=paste(opt$outdir,"/","Homozygote_Model","/",FileName,sep = ""),height=9.5,width=14)

meta1<-metabin(event.e,n.e,event.c,n.c,data=Main,
	method=opt$method,sm="OR",byvar=Main$Country,
	common=TRUE,random=TRUE,
	studlab=paste(Main$Author,Main$Year))

forest(meta1,lab.e="Exposure",lab.c="Control")

dev.off()

# Homozygote model subgroup_analysis_Sample_age forest

event.e<-Sub$GG_case
n.e<-Sub$GG_case+Sub$AA_case
event.c<-Sub$GG_control
n.c<-Sub$GG_control+Sub$AA_control

FileName<-paste("Homozygote model subgroup_analysis_Sample_age.PDF")
pdf(file=paste(opt$outdir,"/","Homozygote_Model","/",FileName,sep = ""),height=9.5,width=14)

meta1<-metabin(event.e,n.e,event.c,n.c,data=Sub,
	method=opt$method,sm="OR",byvar=Sub$Sample_age,
	common=TRUE,random=TRUE,
	studlab=paste(Sub$Author,Sub$Year))

forest(meta1,lab.e="Exposure",lab.c="Control")

dev.off()

# Homozygote model forest_severity

event.e<-Severity$GG_not_severe
n.e<-Severity$GG_not_severe+Severity$AA_not_severe
event.c<-Severity$GG_severe
n.c<-Severity$GG_severe+Severity$AA_severe

meta1<-metabin(event.e,n.e,event.c,n.c,data=Severity,
	method=opt$method,sm="OR", #method = "MH"：指定使用 Mantel-Haenszel法 计算合并效应量
										#特点：适用于研究间效应量同质性较好的情况
										#优势：在小样本或零单元格（zero-cell）情况下表现稳健
										#其他可选方法： 1."Inverse"：逆方差法（默认） 2."Peto"：Peto法，适用于罕见事件（OR接近1） 3."GLMM"：广义线性混合模型
	common=TRUE,random=TRUE,
	studlab=paste(Severity$Author,Severity$Year))

FileName<-paste("Homozygote model forest_severity.PDF")
pdf(file=paste(opt$outdir,"/","Homozygote_Model","/",FileName,sep = ""),height=4.5,width=12)

forest(meta1,lab.e="Not_Severe",lab.c="Severe")

dev.off()

# Homozygote model sensetivity_analysis_severity

FileName<-paste("Homozygote model sensetivity_analysis_",opt$pool,"_severity.PDF")
pdf(file=paste(opt$outdir,"/","Homozygote_Model","/",FileName,sep = ""),height=3.5,width=10)

forest(metainf(meta1,pooled=opt$pool),lab.e="Not_Severe",lab.c="Severe")

dev.off()

# Homozygote model funnel_severity 

FileName<-paste("Homozygote model funnel_severity.PDF")
pdf(file=paste(opt$outdir,"/","Homozygote_Model","/",FileName,sep = ""),height=8,width=8)

funnel(meta1)

dev.off()

# Homozygote model galbraith_severity

FileName<-paste("Homozygote model galbraith_severity.PDF")
pdf(file=paste(opt$outdir,"/","Homozygote_Model","/",FileName,sep = ""),height=8,width=8)

# radial(meta1,level=0.95) #设置置信区间为95%
radial(meta1,text=Severity$Author,level=0.95)

dev.off()


## Recessive model forest 隐性模型

event.e<-Main$GG_case
n.e<-Main$GG_case+Main$AA_AG_case
event.c<-Main$GG_control
n.c<-Main$GG_control+Main$AA_AG_control

meta1<-metabin(event.e,n.e,event.c,n.c,data=Main,
	method=opt$method,sm="OR", #method = "MH"：指定使用 Mantel-Haenszel法 计算合并效应量
										#特点：适用于研究间效应量同质性较好的情况
										#优势：在小样本或零单元格（zero-cell）情况下表现稳健
										#其他可选方法： 1."Inverse"：逆方差法（默认） 2."Peto"：Peto法，适用于罕见事件（OR接近1） 3."GLMM"：广义线性混合模型
	common=TRUE,random=TRUE,
	studlab=paste(Main$Author,Main$Year))

FileName<-paste("Recessive model forest.PDF")
pdf(file=paste(opt$outdir,"/","Recessive_Model","/",FileName,sep = ""),height=4.5,width=12)

forest(meta1,lab.e="Exposure",lab.c="Control")#lab.e:暴露组设置自定义标签（Exposure） lab.c:为对照组设置自定义标签（Control）

dev.off()

# Recessive model sensetivity_analysis

#敏感性分析，注意是用固定效应模型还是随机效应模型
FileName<-paste("Recessive model sensetivity_analysis_",opt$pool,".PDF")
pdf(file=paste(opt$outdir,"/","Recessive_Model","/",FileName,sep = ""),height=3.5,width=10)

forest(metainf(meta1,pooled=opt$pool),lab.e="Exposure",lab.c="Control")

dev.off()

# Recessive model funnel 隐性模型漏斗图

FileName<-paste("Recessive model funnel.PDF")
pdf(file=paste(opt$outdir,"/","Recessive_Model","/",FileName,sep = ""),height=8,width=8)

funnel(meta1)

dev.off()

# Recessive model galbraith 隐性模型甘氏图(定位离群值)

FileName<-paste("Recessive model galbraith.PDF")
pdf(file=paste(opt$outdir,"/","Recessive_Model","/",FileName,sep = ""),height=8,width=8)

# radial(meta1,level=0.95) #设置置信区间为95%
radial(meta1,text=Main$Author,level=0.95)

dev.off()

# Recessive model subgroup_analysis_Country forest

FileName<-paste("Recessive model subgroup_analysis_Country.PDF")
pdf(file=paste(opt$outdir,"/","Recessive_Model","/",FileName,sep = ""),height=9.5,width=14)

meta1<-metabin(event.e,n.e,event.c,n.c,data=Main,
	method=opt$method,sm="OR",byvar=Main$Country,
	common=TRUE,random=TRUE,
	studlab=paste(Main$Author,Main$Year))

forest(meta1,lab.e="Exposure",lab.c="Control")

dev.off()

# Recessive model subgroup_analysis_Sample_age forest

event.e<-Sub$GG_case
n.e<-Sub$GG_case+Sub$AA_AG_case
event.c<-Sub$GG_control
n.c<-Sub$GG_control+Sub$AA_AG_control

FileName<-paste("Recessive model subgroup_analysis_Sample_age.PDF")
pdf(file=paste(opt$outdir,"/","Recessive_Model","/",FileName,sep = ""),height=9.5,width=14)

meta1<-metabin(event.e,n.e,event.c,n.c,data=Sub,
	method=opt$method,sm="OR",byvar=Sub$Sample_age,
	common=TRUE,random=TRUE,
	studlab=paste(Sub$Author,Sub$Year))

forest(meta1,lab.e="Exposure",lab.c="Control")

dev.off()

# Recessive model forest_severity

event.e<-Severity$GG_not_severe
n.e<-Severity$GG_not_severe+Severity$AA_AG_not_severe
event.c<-Severity$GG_severe
n.c<-Severity$GG_severe+Severity$AA_AG_severe

meta1<-metabin(event.e,n.e,event.c,n.c,data=Severity,
	method=opt$method,sm="OR", #method = "MH"：指定使用 Mantel-Haenszel法 计算合并效应量
										#特点：适用于研究间效应量同质性较好的情况
										#优势：在小样本或零单元格（zero-cell）情况下表现稳健
										#其他可选方法： 1."Inverse"：逆方差法（默认） 2."Peto"：Peto法，适用于罕见事件（OR接近1） 3."GLMM"：广义线性混合模型
	common=TRUE,random=TRUE,
	studlab=paste(Severity$Author,Severity$Year))

FileName<-paste("Recessive model forest_severity.PDF")
pdf(file=paste(opt$outdir,"/","Recessive_Model","/",FileName,sep = ""),height=4.5,width=12)

forest(meta1,lab.e="Not_Severe",lab.c="Severe")

dev.off()

# Recessive model sensetivity_analysis_severity

FileName<-paste("Recessive model sensetivity_analysis_",opt$pool,"_severity.PDF")
pdf(file=paste(opt$outdir,"/","Recessive_Model","/",FileName,sep = ""),height=3.5,width=10)

forest(metainf(meta1,pooled=opt$pool),lab.e="Not_Severe",lab.c="Severe")

dev.off()

# Recessive model funnel_severity 

FileName<-paste("Recessive model funnel_severity.PDF")
pdf(file=paste(opt$outdir,"/","Recessive_Model","/",FileName,sep = ""),height=8,width=8)

funnel(meta1)

dev.off()

# Recessive model galbraith_severity

FileName<-paste("Recessive model galbraith_severity.PDF")
pdf(file=paste(opt$outdir,"/","Recessive_Model","/",FileName,sep = ""),height=8,width=8)

# radial(meta1,level=0.95) #设置置信区间为95%
radial(meta1,text=Severity$Author,level=0.95)

dev.off()


##Allele model 等位基因模型

Main<-read.csv(paste(opt$maindata,"allele.csv",sep="_"),header=TRUE,dec=".")
Sub<-read.csv(paste(opt$subdata,"allele.csv",sep="_"),header=TRUE,dec=".")
Severity<-read.csv(paste(opt$severitydata,"allele.csv",sep="_"),header=TRUE,dec=".")

event.e<-Main$G_case
n.e<-Main$Sample_case
event.c<-Main$A_control
n.c<-Main$Sample_control

meta1<-metabin(event.e,n.e,event.c,n.c,data=Main,
	method=opt$method,sm="OR", #method = "MH"：指定使用 Mantel-Haenszel法 计算合并效应量
										#特点：适用于研究间效应量同质性较好的情况
										#优势：在小样本或零单元格（zero-cell）情况下表现稳健
										#其他可选方法： 1."Inverse"：逆方差法（默认） 2."Peto"：Peto法，适用于罕见事件（OR接近1） 3."GLMM"：广义线性混合模型
	common=TRUE,random=TRUE,
	studlab=paste(Main$Author,Main$Year))

FileName<-paste("Allele model forest.PDF")
pdf(file=paste(opt$outdir,"/","Allele_Model","/",FileName,sep = ""),height=4.5,width=12)

forest(meta1,lab.e="Exposure",lab.c="Control")#lab.e:暴露组设置自定义标签（Exposure） lab.c:为对照组设置自定义标签（Control）

dev.off()

# Allele model sensetivity_analysis

#敏感性分析，注意是用固定效应模型还是随机效应模型
FileName<-paste("Allele model sensetivity_analysis_",opt$pool,".PDF")
pdf(file=paste(opt$outdir,"/","Allele_Model","/",FileName,sep = ""),height=3.5,width=10)

forest(metainf(meta1,pooled=opt$pool),lab.e="Exposure",lab.c="Control")

dev.off()

# Allele model funnel 等位基因模型漏斗图

FileName<-paste("Allele model funnel.PDF")
pdf(file=paste(opt$outdir,"/","Allele_Model","/",FileName,sep = ""),height=8,width=8)

funnel(meta1)

dev.off()

# Allele model galbraith 等位基因模型甘氏图(定位离群值)

FileName<-paste("Allele model galbraith.PDF")
pdf(file=paste(opt$outdir,"/","Allele_Model","/",FileName,sep = ""),height=8,width=8)

# radial(meta1,level=0.95) #设置置信区间为95%
radial(meta1,text=Main$Author,level=0.95)

dev.off()

# Allele model subgroup_analysis_Country forest

FileName<-paste("Allele model subgroup_analysis_Country.PDF")
pdf(file=paste(opt$outdir,"/","Allele_Model","/",FileName,sep = ""),height=9.5,width=14)

meta1<-metabin(event.e,n.e,event.c,n.c,data=Main,
	method=opt$method,sm="OR",byvar=Main$Country,
	common=TRUE,random=TRUE,
	studlab=paste(Main$Author,Main$Year))

forest(meta1,lab.e="Exposure",lab.c="Control")

dev.off()

# Allele model subgroup_analysis_Sample_age forest

event.e<-Sub$G_case
n.e<-Sub$Sample_case
event.c<-Sub$A_control
n.c<-Sub$Sample_control

FileName<-paste("Allele model subgroup_analysis_Sample_age.PDF")
pdf(file=paste(opt$outdir,"/","Allele_Model","/",FileName,sep = ""),height=9.5,width=14)

meta1<-metabin(event.e,n.e,event.c,n.c,data=Sub,
	method=opt$method,sm="OR",byvar=Sub$Sample_age,
	common=TRUE,random=TRUE,
	studlab=paste(Sub$Author,Sub$Year))

forest(meta1,lab.e="Exposure",lab.c="Control")

dev.off()

# Allele model forest_severity

event.e<-Severity$G_not_severe
n.e<-Severity$Sample_not_severe
event.c<-Severity$A_severe
n.c<-Severity$Sample_severe

meta1<-metabin(event.e,n.e,event.c,n.c,data=Severity,
	method=opt$method,sm="OR", #method = "MH"：指定使用 Mantel-Haenszel法 计算合并效应量
										#特点：适用于研究间效应量同质性较好的情况
										#优势：在小样本或零单元格（zero-cell）情况下表现稳健
										#其他可选方法： 1."Inverse"：逆方差法（默认） 2."Peto"：Peto法，适用于罕见事件（OR接近1） 3."GLMM"：广义线性混合模型
	common=TRUE,random=TRUE,
	studlab=paste(Severity$Author,Severity$Year))

FileName<-paste("Allele model forest_severity.PDF")
pdf(file=paste(opt$outdir,"/","Allele_Model","/",FileName,sep = ""),height=4.5,width=12)

forest(meta1,lab.e="Not_Severe",lab.c="Severe")

dev.off()

# Allele model sensetivity_analysis_severity

FileName<-paste("Allele model sensetivity_analysis_",opt$pool,"_severity.PDF")
pdf(file=paste(opt$outdir,"/","Allele_Model","/",FileName,sep = ""),height=3.5,width=10)

forest(metainf(meta1,pooled=opt$pool),lab.e="Not_Severe",lab.c="Severe")

dev.off()

# Allele model funnel_severity 

FileName<-paste("Allele model funnel_severity.PDF")
pdf(file=paste(opt$outdir,"/","Allele_Model","/",FileName,sep = ""),height=8,width=8)

funnel(meta1)

dev.off()

# Allele model galbraith_severity

FileName<-paste("Allele model galbraith_severity.PDF")
pdf(file=paste(opt$outdir,"/","Allele_Model","/",FileName,sep = ""),height=8,width=8)

# radial(meta1,level=0.95) #设置置信区间为95%
radial(meta1,text=Severity$Author,level=0.95)

dev.off()

