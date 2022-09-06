# file_snp = 'test.pileup'
# file_blacklist = 'blacklist.final.bed'


library(tidyverse)
library(data.table)
library(GenomicRanges)

args = commandArgs(trailingOnly = T)
file_snp = args[1]
file_blacklist = args[2]

out.file_snp_rmblk = str_replace(file_snp,'.gz','.no_blk.gz')


# read blacklist file
df_blk = read.table(file_blacklist, header = F)
colnames(df_blk) = c('chr','start','end','name')
g_blk = makeGRangesFromDataFrame(df_blk)


# read count file
df_snpcount = fread(file_snp, sep = ',', header = T)
df_snpcount$start = df_snpcount$Position - 1
df_snpcount$end = df_snpcount$Position
if (grepl(pattern = 'chr|Chr', x = df_snpcount[1,1]) == TRUE){
	df_snpcount$Chromosome = str_replace(string = df_snpcount$Chromosome, pattern = 'chr|Chr',replace = '')
}

#########3
# make granges
g_snpcount = makeGRangesFromDataFrame(df_snpcount, keep.extra.columns = T)
# filter snp file
g_snpcount_new = g_snpcount[-queryHits(findOverlaps(g_snpcount, g_blk, type = 'any')),]
df_snpcount_new = data.frame(g_snpcount_new)
df_snpcount_new$start = NULL
df_snpcount_new$end = NULL
df_snpcount_new$width = NULL
df_snpcount_new$strand = NULL
df_snpcount_new = df_snpcount_new %>% dplyr::rename('Chromosome' = 'seqnames')


# write file to the output
write.csv(df_snpcount_new, quote = F, row.names = F , file = gzfile(out.file_snp_rmblk))

