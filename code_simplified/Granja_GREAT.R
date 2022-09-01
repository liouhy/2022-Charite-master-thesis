# perform GO enrichment
.libPaths('/fast/home/h/hliou/miniconda3/envs/r-kernel/lib/R/library')
library(rGREAT)

setwd('/fast/AG_Ohler/hyliou/mtproject/code')

bg = read.table('../processed/Granja_ATAC_bed.bed', sep='\t')
bg = bg[grep('chr', bg[,1]),]

great_go = function(acc='acc') {
  
  dir = paste0('../processed/Granja/',acc,'_regions/')
  fl = list.files(dir)
  for (file in fl) {
    
    err = FALSE
    tryCatch({test_acc = read.table(paste0(dir,file), sep='\t')},
             error=function(e){err<<-TRUE})
    if (err){next}
    
    if (nrow(test_acc)==0) {next}
    
    test_acc = test_acc[grep('chr', test_acc[,1]),]
    
    job_acc = submitGreatJob(test_acc, bg, species='hg19')
    
    tb_acc = getEnrichmentTables(job_acc, ontology=c('GO Biological Process'))
    
    df_acc = tb_acc[[1]]
    df_acc = df_acc[order(df_acc['Hyper_Adjp_BH']),]
    write.table(df_acc,paste0('../Granja/GREAT/',acc,strsplit(file,"[.]")[[1]][1],'.tsv'),sep='\t')
  }
  
}
# GREAT GO enrichment for enriched regions of each reference cell
great_go('acc')

great_go('inacc')
