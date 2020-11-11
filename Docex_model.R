OAK_Docex = read.csv("OAK_Docex.csv",header = T)
OAK_Docex = OAK_Docex[,c("PtID","bTMB_total","bTMB_sense","UMS_OS","PFS",
                                       "PFS.CNSR","OS","OS.CNSR","SEX","HIST","TOBHX",
                                       "METSITES","KRASMUT","EGFRMUT","EML4MUT")]
colnames(OAK_Docex)[4] = "UMS"
OAK_Docex$DRIVER = apply(OAK_Docex[,c("KRASMUT","EGFRMUT","EML4MUT")], 1, function(x){
  status = "NEGATIVE"
  if(sum(x == "POSITIVE") >=1){status = "POSITIVE"}
  return(status)
})


library(survival)
library(survminer)
library(ggplot2)
library(ggsci)
library(reshape2)
find_threshold = function(df, indicator = "bTMB_total"){
  require(survival, quietly = T)
  pvalue = data.frame()
  for(i in 1:30){
    df$status = "Low"
    df$status[df[,indicator] > i] = "High"
    df$status = factor(df$status, levels = c("Low","High"))
    if(nlevels(df$status) < 2){next}
    coxph_res = summary(coxph(Surv(PFS, PFS.CNSR) ~ status, data = df))
    PFS_pvalue = coxph_res$coefficients[5]
    coxph_res = summary(coxph(Surv(OS, OS.CNSR) ~ status, data = df))
    OS_pvalue = coxph_res$coefficients[5]
    hr = coxph_res$coefficients[2]
    pvalue = rbind(pvalue, data.frame(Threshold = i, min_Sample = min(table(df$status)),hazard_ratio = hr,
                                      PFS = PFS_pvalue, OS = OS_pvalue, stringsAsFactors = F))
  }
  th = pvalue[pvalue$min_Sample >=nrow(df)/4,]
  ind = which(th$OS == min(th$OS))[1]
  return(list(Pvalues = pvalue, optimal_th = th$Th[ind], hazard_ratio = th$hazard_ratio[ind]))
}



## total bTMB
bTMB_pvalue = find_threshold(OAK_Docex, indicator = "bTMB_total")
ggplot(melt(bTMB_pvalue$Pvalues[,c("Threshold", "PFS", "OS")], id = "Threshold"), mapping = aes(x=Threshold, y=value, color=variable)) + 
  geom_line(size = 1) + geom_point(size = 2) + scale_y_log10() + ylab("p value") +
  scale_color_npg() + geom_vline(xintercept = bTMB_pvalue$optimal_th, linetype = 2) +
  theme_bw() + theme(axis.title = element_text(size=18, color = "black"),
                     axis.text = element_text(size=16, color = "black"),
                     legend.text = element_text(size=16),
                     legend.title = element_blank())
## sense bTMB
sbTMB_pvalue = find_threshold(OAK_Docex, indicator = "bTMB_sense")
ggplot(melt(sbTMB_pvalue$Pvalues[,c("Threshold", "PFS", "OS")], id = "Threshold"), mapping = aes(x=Threshold, y=value, color=variable)) + 
  geom_line(size = 1) + geom_point(size = 2) + scale_y_log10() + ylab("p value") +
  scale_color_npg() + geom_vline(xintercept = sbTMB_pvalue$optimal_th, linetype = 2) +
  theme_bw() + theme(axis.title = element_text(size=18, color = "black"),
                     axis.text = element_text(size=16, color = "black"),
                     legend.text = element_text(size=16),
                     legend.title = element_blank())

## UMs
UMS_pvalue = find_threshold(OAK_Docex, indicator = "UMS")
ggplot(melt(UMS_pvalue$Pvalues[,c("Threshold", "PFS", "OS")], id = "Threshold"), mapping = aes(x=Threshold, y=value, color=variable)) + 
  geom_line(size = 1) + geom_point(size = 2) + scale_y_log10() + ylab("p value") +
  scale_color_npg() + geom_vline(xintercept = UMS_pvalue$optimal_th, linetype = 2) +
  theme_bw() + theme(axis.title = element_text(size=18, color = "black"),
                     axis.text = element_text(size=16, color = "black"),
                     legend.text = element_text(size=16),
                     legend.title = element_blank())


hazard_ratio = data.frame()
optimal_threshold = data.frame()
for(f in c("SEX","HIST","TOBHX","METSITES","DRIVER")){
  groups = unique(OAK_Docex[, f])
  hr = data.frame()
  th = data.frame()
  for(g in groups){
    df = OAK_Docex[OAK_Docex[,f] == g,]
    sub_bTMB_pvalue = find_threshold(df, indicator = "bTMB_total")
    sub_sbTMB_pvalue = find_threshold(df, indicator = "bTMB_sense")
    sub_UMS_pvalue = find_threshold(df, indicator = "UMS")
    
    hr =rbind(hr, data.frame(bTMB = sub_bTMB_pvalue$hazard_ratio, sbTMB = sub_sbTMB_pvalue$hazard_ratio,
                             UMS = sub_UMS_pvalue$hazard_ratio, stringsAsFactors = F)) 
    th =rbind(th, data.frame(bTMB = sub_bTMB_pvalue$optimal_th, sbTMB = sub_sbTMB_pvalue$optimal_th,
                             UMS = sub_UMS_pvalue$optimal_th, stringsAsFactors = F))
  }
  hr$group = groups
  hr$factors = f
  th$group = groups
  th$factors = f
  hazard_ratio = rbind(hazard_ratio, hr)
  optimal_threshold = rbind(optimal_threshold, th)
}

hazard_ratio = rbind(hazard_ratio, data.frame(bTMB = bTMB_pvalue$hazard_ratio, sbTMB = sbTMB_pvalue$hazard_ratio,
                                              UMS = UMS_pvalue$hazard_ratio, group = "total", factors = "total",
                                              stringsAsFactors = F))
optimal_threshold = rbind(optimal_threshold, data.frame(bTMB = bTMB_pvalue$optimal_th, sbTMB = sbTMB_pvalue$optimal_th,
                                                        UMS = UMS_pvalue$optimal_th, group = "total", factors = "total",
                                                        stringsAsFactors = F))

score_patients = c()
for(i in 1:nrow(OAK_Docex)){
  SEX = OAK_Docex$SEX[i]
  HIST = OAK_Docex$HIST[i]
  TOBHX = OAK_Docex$TOBHX[i]
  METSITES = OAK_Docex$METSITES[i]
  DRIVER = OAK_Docex$DRIVER[i]
  index = c(SEX, TOBHX, HIST, METSITES, DRIVER)
  # bTMB
  score1 = sum(as.numeric(OAK_Docex$bTMB_total[i] > 
                            optimal_threshold$bTMB[optimal_threshold$group %in% index]) * 
                 sum(hazard_ratio$bTMB[hazard_ratio$group %in% index]))
  
  score2 = sum(as.numeric(OAK_Docex$bTMB_sense[i] > 
                            optimal_threshold$sbTMB[optimal_threshold$group %in% index]) * 
                 sum(hazard_ratio$sbTMB[hazard_ratio$group %in% index]))
  
  score3 = sum(as.numeric(OAK_Docex$UMS[i] > 
                            optimal_threshold$UMS[optimal_threshold$group %in% index]) * 
                 sum(hazard_ratio$UMS[hazard_ratio$group %in% index]))
  
  score_patients = c(score_patients, sum(c(score1, score2, score3)))
}
OAK_Docex$score = score_patients


OAK_Docex$status = "LOW"
OAK_Docex$status[OAK_Docex$score>=median(OAK_Docex$score)] = "HIGH"
OAK_Docex$status = factor(OAK_Docex$status, levels = c("LOW","HIGH"))

cox_res = coxph(Surv(OS, OS.CNSR) ~ status, data = OAK_Docex)
fit = survfit(Surv(OS, OS.CNSR) ~ status, data = OAK_Docex)
ggsurvplot(fit, data = OAK_Docex,
           conf.int = TRUE,#添加置信区间
           pval = TRUE,#添加P值
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           linetype = "strata",
           palette = c("#E7B800","#2E9FDF"),
           legend = "bottom",
           legend.title = "RISK",
           legend.labs = c("LOW","HIGH"))

write.csv(OAK_Docex, file = "OAK_Docex_res.csv", row.names = F, quote = F)

#### validation
POPLAR_Docex = read.csv("POPLAR_Docex.csv",header = T)
POPLAR_Docex = POPLAR_Docex[,c("PtID","bTMB_total","bTMB_sense","UMS_OS","PFS",
                                             "PFS.CNSR","OS","OS.CNSR","SEX","HIST","TOBHX",
                                             "METSITES","KRASMUT","EGFRMUT","EML4MUT")]
colnames(POPLAR_Docex)[4] = "UMS"
POPLAR_Docex$DRIVER = apply(POPLAR_Docex[,c("KRASMUT","EGFRMUT","EML4MUT")], 1, function(x){
  status = "NEGATIVE"
  if(sum(x == "POSITIVE") >=1){status = "POSITIVE"}
  return(status)
})

score_patients = c()
for(i in 1:nrow(POPLAR_Docex)){
  SEX = POPLAR_Docex$SEX[i]
  HIST = POPLAR_Docex$HIST[i]
  TOBHX = POPLAR_Docex$TOBHX[i]
  METSITES = POPLAR_Docex$METSITES[i]
  DRIVER = POPLAR_Docex$DRIVER[i]
  index = c(SEX, TOBHX, HIST, METSITES, DRIVER)
  # bTMB
  
  score1 = sum(as.numeric(POPLAR_Docex$bTMB_total[i] > 
                            optimal_threshold$bTMB[optimal_threshold$group %in% index]) * 
                 sum(hazard_ratio$bTMB[hazard_ratio$group %in% index]))
  
  score2 = sum(as.numeric(POPLAR_Docex$bTMB_sense[i] > 
                            optimal_threshold$sbTMB[optimal_threshold$group %in% index]) * 
                 sum(hazard_ratio$sbTMB[hazard_ratio$group %in% index]))
  
  score3 = sum(as.numeric(POPLAR_Docex$UMS[i] > 
                            optimal_threshold$UMS[optimal_threshold$group %in% index]) * 
                 sum(hazard_ratio$UMS[hazard_ratio$group %in% index]))
  
  score_patients = c(score_patients, sum(c(score1, score2, score3)))
}
POPLAR_Docex$score = score_patients

POPLAR_Docex$status = "LOW"
POPLAR_Docex$status[POPLAR_Docex$score >= median(OAK_Docex$score)] = "HIGH"
POPLAR_Docex$status = factor(POPLAR_Docex$status, levels = c("LOW","HIGH"))

cox_res = coxph(Surv(OS, OS.CNSR) ~ status, data = POPLAR_Docex)
fit = survfit(Surv(OS, OS.CNSR) ~ status, data = POPLAR_Docex)
ggsurvplot(fit, data = POPLAR_Docex,
           conf.int = TRUE,#添加置信区间
           pval = TRUE,#添加P值
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           linetype = "strata",
           palette = c("#E7B800","#2E9FDF"),
           legend = "bottom",
           legend.title = "RISK",
           legend.labs = c("LOW","HIGH"))

write.csv(POPLAR_Docex, file = "POPLAR_Docex_res.csv", row.names = F, quote = F)

