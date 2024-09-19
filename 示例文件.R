library(SMRinR)

#在线获取工具变量
UKB_PPP_pQTLs <- get_cis_pQTL_from_UKB_PPP()
DECODE_pQTLs <- get_cis_pQTL_from_DECODE()
INTERVAL_pQTLs <- get_cis_pQTL_from_INTERVAL()
Fenland_pQTLs <- get_cis_pQTL_from_Fenland()
FinnGen_Olink_pQTLs <- get_cis_pQTL_from_FinnGen_Olink()
FinnGen_Somascan_pQTLs <- get_cis_pQTL_from_FinnGen_Somascan()
# 注意！！ exposure列为蛋白ID（唯一识别码），info.exposure为蛋白名称，gene.exposure为蛋白所对应的基因名称
# 注意！！ exposure列为蛋白ID（唯一识别码），info.exposure为蛋白名称，gene.exposure为蛋白所对应的基因名称

# 工具变量筛选原则
show_details()


# 在线获取所有pQTL的结局信息
# 按需要修改 UKB_PPP_pQTLs 这个变量
outcomes <- TwoSampleMR::extract_outcome_data(UKB_PPP_pQTLs$SNP, 'ukb-b-19953', proxies = F)
dat <- TwoSampleMR::harmonise_data(UKB_PPP_pQTLs, outcomes)

# 本地结局（以finn为例）
all_outcome <- data.table::fread('G:/迅雷下载/finngen_R9_J10_ARDS.gz')
outcomes <- all_outcome[all_outcome$rsids %in% IVs$SNP, ]
outcomes <- as.data.frame(outcomes)
outcomes <-  TwoSampleMR::format_data(
  outcomes,
  type = "outcome",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  eaf_col = "af_alt",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  pval_col = "pval"
)
dat <- TwoSampleMR::harmonise_data(IVs, outcomes)



# SMR
result <- SMR(dat)
result <- result[result$mr_keep == TRUE, ]
result$pval_adj_SMR <- p.adjust(result$pval_SMR, method = 'fdr')


# SMR结果可视化
data <- result
data$log_pval <- -log10(data$pval_adj_SMR)
data$significance <- with(data, ifelse(pval_adj_SMR < 0.05 & beta_SMR > 0, "Positively Associated",
                                       ifelse(pval_adj_SMR < 0.05 & beta_SMR < 0, "Negatively Associated", "Not Significant")))
color_values <- c("Positively Associated" = "red", "Negatively Associated" = "blue", "Not Significant" = "grey")


library(ggplot2)
library(ggrepel)
ggplot(data, aes(x = beta_SMR, y = log_pval)) +
  geom_point(aes(color = significance), alpha = 0.8, size = 3) +
  scale_color_manual(values = color_values) +
  theme_minimal() +
  labs(
    title = "Volcano Plot of SMR Results",
    x = "Beta (Effect Size)",
    y = "-log10(Adjusted P-Value)",
    color = "Association"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +  # 
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 1) +  # 
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 10, hjust = 0.5),
    panel.grid = element_blank()
  ) +
  geom_text_repel(data = subset(data, significance != "Not Significant"),
                  aes(label = exposure),
                  size = 3.5,
                  box.padding = 0.3,
                  point.padding = 0.2,
                  max.overlaps = 10)





# HEIDI test与共定位（一次性做完两个检测）
# 注意：无需对所有的SMR结果都进行HEIDI test与共定位，对多重校正之后p<0.05的结果进行HEIDI test与共定位，节约时间。
# 对所选的top_SNP上下1 Mb的SNP截取下来，并提取结局信息，以上步骤均使用TwosampleMR即可。
# 现在支持对deCODE, UKB-PPP, Fenland, FinnGen_Olink, FinnGen_Somascan的在线获取,INTERvAL在opengwas上有

# !!!
# 需要传入id，id为大家填写表单时候的邮箱
# 需设置ieu token，heidi test会用到

# deCODE示例
# protein_name为DECODE_pQTLs的exposure列的名字
# ！！exposure列
# ！！exposure列
# start和end请大家自行修改
protein_name <- "GFAP"
chrom_target <- 10
pos_start <- 100000
pos_end <- 302000
# 需要传入id，id为大家填写表单时候的邮箱
id <- 'XXX'
IVs <- get_deCODE_data(protein_name, chrom_target, pos_start, pos_end,id)
IVs <- TwoSampleMR::format_data(
  IVs,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  snp_col = "rsids",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col = "effectAllele",
  other_allele_col = "otherAllele",
  pval_col = "Pval",
  samplesize_col = "N",
  chr_col = "Chrom",
  pos_col = "Pos",
  log_pval = FALSE
)
outcomes <- TwoSampleMR::extract_outcome_data(IVs$SNP, 'ukb-b-19953', proxies = F)
dat <- TwoSampleMR::harmonise_data(IVs, outcomes)

library(coloc)
# 共定位
D1<-list(type="quant",
         beta=dat$beta.exposure,
         varbeta=dat$se.exposure^2,
         N=838,
         snp=dat$SNP,
         sdY=1)

D2<-list(type="quant",
         beta=dat$beta.outcome,
         varbeta=dat$se.outcome^2,
         N=838,
         snp=dat$SNP,
         sdY=1)

coloc<-coloc.abf( D1, D2, p1 = 1e-04,p2 = 1e-04,p12 = 1e-05)

# HEIDI test
top_SNP <- 'rs141533188'
HEIDI_test(dat,top_SNP)




# UKB-PPP示例
# protein_name为UKB_PPP_pQTLs的exposure列的名字
# ！！exposure列
# ！！exposure列
# start和end请大家自行修改
library(magrittr)
protein_name <- "A1BG"
chrom_target <- 11
pos_start <- 100000
pos_end <- 602000
# 需要传入id，id为大家填写表单时候的邮箱
id <- 'XXX'
IVs <- get_UKB_PPP_data(protein_name, chrom_target, pos_start, pos_end,id)
# 这里需要转rsid,构建一个df，注意 第一列是ID，第二列是CHROM，顺序不能调换。
df <- IVs[,c(3,1)]
mapping_result <- rsid_mapping(df,id)
IVs_with_rsid <- merge(IVs, mapping_result,by = 'ID')
# 整理格式
IVs <- TwoSampleMR::format_data(
  IVs_with_rsid,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  snp_col = "rsid",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "LOG10P",
  samplesize_col = "N",
  chr_col = "CHROM",
  pos_col = "GENPOS",
  log_pval = TRUE
)
outcomes <- TwoSampleMR::extract_outcome_data(IVs$SNP, 'ukb-b-19953', proxies = F)
dat <- TwoSampleMR::harmonise_data(IVs, outcomes)
# 共定位和HEIDI test代码与之前一致


# Fenland示例
# ！！！！！seqid为Fenland_pQTLs的exposure列的名字
# ！！exposure列
# ！！exposure列
# start和end请大家自行修改
seqid <- "SeqId_9191_8"
chromosome <- 20
pos_start <- 100000
pos_end <- 302000
# 需要传入id，id为大家填写表单时候的邮箱
id <- 'XXX'
IVs <- get_Fenland_data(seqid, chromosome, pos_start, pos_end, id)
IVs <- TwoSampleMR::format_data(
  IVs,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  snp_col = "rsid",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = 'Freq1',
  pval_col = "Pvalue",
  samplesize_col = "TotalSampleSize",
  chr_col = "chr",
  pos_col = "pos",
  log_pval = FALSE
)
outcomes <- TwoSampleMR::extract_outcome_data(IVs$SNP, 'ukb-b-19953', proxies = F)
dat <- TwoSampleMR::harmonise_data(IVs, outcomes)
# 共定位和HEIDI test代码与之前一致


# FinnGen_Olink示例
# ！！！！！seqid为FinnGen_Olink_pQTLs的exposure列的名字
# ！！exposure列
# ！！exposure列
# start和end请大家自行修改
source <- "Olink"
target <- "AAMDC"
chr_value <- 5
pos_min <- 1000000
pos_max <- 1200000
# 需要传入id，id为大家填写表单时候的邮箱
id <- 'XXX'
# 调用函数并获取结果
IVs <- get_FinnGen_data(source, target, chr_value, pos_min, pos_max)
# 这里需要转rsid,构建一个df，注意 第一列是ID，第二列是CHROM，顺序不能调换。
df <- IVs[,c(4,5)]
mapping_result <- rsid_mapping2(df,id)
mapping_result$variant <- paste0("chr", mapping_result$variant)
mapping_result$variant <- gsub(":", "_", mapping_result$variant)
IVs_with_rsid <- merge(IVs, mapping_result,by.x = 'ID',by.y = 'variant')
# 整理格式
IVs <- TwoSampleMR::format_data(
  IVs_with_rsid,
  type = "exposure",
  header = TRUE,
  snp_col = "rsid",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "ALT_FREQ",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P",
  samplesize_col = "N",
  chr_col = "CHR",
  pos_col = "POS",
)
outcomes <- TwoSampleMR::extract_outcome_data(IVs$SNP, 'ukb-b-19953', proxies = F)
dat <- TwoSampleMR::harmonise_data(IVs, outcomes)
# 共定位和HEIDI test代码与之前一致


# FinnGen_Somascan示例
# ！！！！！seqid为FinnGen_Somascan_pQTLs的exposure列的名字
# ！！exposure列
# ！！exposure列
# start和end请大家自行修改
source <- "Somascan"
target <- "seq.10037.98"
chr_value <- 5
pos_min <- 1000000
pos_max <- 1200000
# 需要传入id，id为大家填写表单时候的邮箱
id <- 'XXX'
# 调用函数并获取结果
IVs <- get_FinnGen_data(source, target, chr_value, pos_min, pos_max)
# 这里需要转rsid,构建一个df，注意 第一列是ID，第二列是CHROM，顺序不能调换。
df <- IVs[,c(4,5)]
mapping_result <- rsid_mapping2(df,id)
mapping_result$variant <- paste0("chr", mapping_result$variant)
mapping_result$variant <- gsub(":", "_", mapping_result$variant)
IVs_with_rsid <- merge(IVs, mapping_result,by.x = 'ID',by.y = 'variant')
# 整理格式
IVs <- TwoSampleMR::format_data(
  IVs_with_rsid,
  type = "exposure",
  header = TRUE,
  snp_col = "rsid",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "ALT_FREQ",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P",
  samplesize_col = "N",
  chr_col = "CHR",
  pos_col = "POS",
)
outcomes <- TwoSampleMR::extract_outcome_data(IVs$SNP, 'ukb-b-19953', proxies = F)
dat <- TwoSampleMR::harmonise_data(IVs, outcomes)
# 共定位和HEIDI test代码与之前一致




# 推荐使用smr.exe完成分析，因为共定位和HEIDI 需要用smr.exe提取原始数据，但是可以用R包先进行初筛

# GTEx eqtl
GTEx_eQTLs <- read.csv("eqtl.csv")
GTEx_eQTLs <- GTEx_eQTLs[GTEx_eQTLs$pval_nominal <= GTEx_eQTLs$pval_nominal_threshold, ]
GTEx_eQTLs <- GTEx_eQTLs[GTEx_eQTLs$chr == GTEx_eQTLs$gene_chr, ]
GTEx_eQTLs <- GTEx_eQTLs[GTEx_eQTLs$variant_pos >=  GTEx_eQTLs$gene_start - 2000000 & GTEx_eQTLs$variant_pos <=  GTEx_eQTLs$gene_end + 2000000, ]
GTEx_eQTLs_IVs <- TwoSampleMR::format_data(
  GTEx_eQTLs,
  type = "exposure",
  phenotype_col = "gene_id",
  snp_col = "variant_id",
  beta_col = "slope",
  se_col = "slope_se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  pval_col = "pval_nominal",
  chr_col = "chr",
  pos_col = "variant_pos",
  log_pval = FALSE
)



# GTEx sqtl
GTEx_sQTLs <- read.csv("sqtl.csv")
GTEx_sQTLs <- GTEx_sQTLs[GTEx_sQTLs$pval_nominal <= GTEx_sQTLs$pval_nominal_threshold, ]
GTEx_sQTLs <- GTEx_sQTLs[GTEx_sQTLs$chr == GTEx_sQTLs$gene_chr, ]
GTEx_sQTLs <- GTEx_sQTLs[GTEx_sQTLs$variant_pos >=  GTEx_sQTLs$gene_start - 2000000 & GTEx_sQTLs$variant_pos <=  GTEx_sQTLs$gene_end + 2000000, ]
GTEx_sQTLs_IVs <- TwoSampleMR::format_data(
  GTEx_sQTLs,
  type = "exposure",
  phenotype_col = "gene_id",
  snp_col = "variant_id",
  beta_col = "slope",
  se_col = "slope_se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  pval_col = "pval_nominal",
  chr_col = "chr",
  pos_col = "variant_pos",
  log_pval = FALSE
)

