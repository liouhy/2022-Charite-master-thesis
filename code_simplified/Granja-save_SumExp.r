# Read the RDS file and save data into other formats
library(SummarizedExperiment)
library(Matrix)

setwd('../')
se = readRDS('raw/scATAC-Healthy-Hematopoiesis-191120.rds')
c_m = assays(se)[[1]]
r_d = rowData(se)
c_d = colData(se)

writeMM(c_m, 'raw/scATAC-Healthy-Hematopoiesis-191120-counts.mtx')
write.csv(r_d, 'raw/scATAC-Healthy-Hematopoiesis-191120-rows.csv', row.names=F)
write.csv(c_d, 'raw/scATAC-Healthy-Hematopoiesis-191120-cols.csv', row.names=F)