library(tidyverse)

args = commandArgs(trailingOnly = T)

file_dbsnp = args[1]
file_XY = args[2]

file_dbsnp_new = str_replace(string = file_dbsnp, pattern = '.gz', replace = '.add_mid_exon.gz')

list_col_types = list(col_character(), col_double(), col_character(), col_character(), col_double(), col_double(), col_double(), col_double(), col_double(), col_double(), col_double(), col_double())

df.dbsnp = read_csv(file_dbsnp, col_types = list_col_types)
df.XY = read_csv(file_XY, col_types = list_col_types)

df.dbsnp$type = 'dbsnp'
df.XY$Ref = '.'
df.XY$Alt = '.'
df.XY$type = 'mid_exon'

# merge them together: (1) select dbsnp if they're duplicated (2)process X & Y separately
df.dbsnp.auto = df.dbsnp %>% filter(!(Chromosome %in% c('X','Y')))
df.dbsnp.XY = df.dbsnp %>% filter((Chromosome %in% c('X','Y')))

df.merge.XY = rbind(df.dbsnp.XY, df.XY) 
df.merge.XY$type = factor(df.merge.XY$type, levels = c('dbsnp','mid_exon'))
df.merge.XY = df.merge.XY %>% arrange(Chromosome, Position, type) %>% distinct(Chromosome, Position, .keep_all = T) 


# rbind: autosome & XY
df.merge = rbind(df.dbsnp.auto , df.merge.XY)
df.merge %>% select(-type) %>% write_csv(file_dbsnp_new)





