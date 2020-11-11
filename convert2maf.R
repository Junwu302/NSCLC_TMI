mutation_df = read.table("OAK_POPLAR_btmb_variants.txt",header = T, sep='\t',stringsAsFactors = F)
mutation_df = mutation_df[,c(1,3,5,6,10,11,12,13,16)]
colnames(mutation_df) = c("Tumor_Sample_Barcode","Chromosome","Reference_Allele","Tumor_Seq_Allele2",
                          "Start_Position","End_Position","Hugo_Symbol","i_transcript_name","Variant_Classification")
mutation_df$Variant_Type = "SNP"
mutation_df$Variant_Classification[mutation_df$Variant_Classification=="missense"] = "Missense_Mutation"
mutation_df$Variant_Classification[mutation_df$Variant_Classification=="nonsense"] = "Nonsense_Mutation"
mutation_df$Variant_Classification[mutation_df$Variant_Classification=="splice"] = "Splice_Site"
mutation_df$Variant_Classification[mutation_df$Variant_Classification=="synonymous"] = "synonymous"

#Clinical information
Docex_ptID = read.table("Docex_oncoplot_PtID_selected.txt", header = T, sep='\t',stringsAsFactors = F)
Atezolizumab_ptID = read.table("Atezolizumab_oncoplot_PtID_selected.txt", header = T, sep='\t',stringsAsFactors = F)
Docex_ptID$group = "Docex"
Atezolizumab_ptID$group = "Atezolizumab"
ptID = rbind(Docex_ptID, Atezolizumab_ptID)

ptID$TOBHX[ptID$TOBHX == 'NEVER'] = "NO"
ptID$TOBHX[ptID$TOBHX %in% c('CURRENT',"PREVIOUS")] = "YES"

ptID$METSITES[Docex_ptID$METSITES %in% c("1","2","< 3")] = "< 3"
ptID$METSITES[!Docex_ptID$METSITES %in% c("1","2","< 3")] = "≥ 3"

ptID$DRIVER = apply(ptID[,6:8],1,function(x){
  y = "NEGATIVE"
  if(sum(x=="POSITIVE")>0){
    y="POSITIVE"
  }
  return(y)
})
ptID = ptID[,-(6:8)]
colnames(ptID)[1] = "Tumor_Sample_Barcode"




#selected patients
Docex_mutation = mutation_df[mutation_df$Tumor_Sample_Barcode %in% ptID$Tumor_Sample_Barcode[ptID$group == "Docex"],]
Atezolizumab_mutation = mutation_df[mutation_df$Tumor_Sample_Barcode %in% ptID$Tumor_Sample_Barcode[ptID$group == "Atezolizumab"],]

write.table(Docex_mutation, file = "Docex_mutation.txt",sep='\t',row.names = F, col.names = T, quote = F)
write.table(Atezolizumab_mutation, file = "Atezolizumab_mutation.txt",sep='\t',row.names = F, col.names = T, quote = F)


Docex_laml = read.maf("Docex_mutation.txt",clinicalData = ptID[ptID$group == "Docex",])
Atezolizumab_laml = read.maf("Atezolizumab_mutation.txt",clinicalData = ptID[ptID$group == "Atezolizumab",])

pdf("Docex_mafsummary.pdf",width = 10, height = 6)
plotmafSummary(maf = Docex_laml, rmOutlier = F, addStat = 'median', 
               dashboard = TRUE, titvRaw = FALSE)
dev.off()

pdf("Docex_oncoplot.pdf", width = 14, height = 5)
oncoplot(maf = Docex_laml, top = 20,
         clinicalFeatures=c("SEX","HIST","TOBHX","METSITES","DRIVER"))
dev.off()

pdf("Atezolizumab_mafsummary.pdf",width = 10, height = 6)
plotmafSummary(maf = Atezolizumab_laml, rmOutlier = F, addStat = 'median', 
               dashboard = TRUE, titvRaw = FALSE)
dev.off()

pdf("Atezolizumab_oncoplot.pdf", width = 14, height = 5)
oncoplot(maf = Atezolizumab_laml, top = 20,
         clinicalFeatures=c("SEX","HIST","TOBHX","METSITES","DRIVER"))
dev.off()

