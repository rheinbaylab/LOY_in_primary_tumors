library(devtools)
library(tidyverse)


library('facetschrY4')

genome_build = "hg19_Y"
load(system.file("ext_genomesize", "genomesize.rda", package = "facetschrY4"))
load(system.file("ext_genomesize", "hg19_Y.rda", package = "facetschrY4"))
df_genomesize = list.genomesize[[genome_build]]



args = commandArgs(trailingOnly=TRUE)

facets_rdata = args[1]
output_prefix = args[2]

# facets_rdata = "test-TP-NB.facets.RData"
# output_prefix = "test-TP-NB"



filename_CNY_stat = paste(output_prefix, '.facets_LOY.stat.tsv', sep='')
filename_CNY_is_LOX = paste(output_prefix, '.tmp.facets_LOY.isLOX.txt', sep='')
filename_CNY_is_LOY = paste(output_prefix, '.tmp.facets_LOY.isLOY.txt', sep='')
filename_CNY_inferred_sex = paste(output_prefix, '.tmp.facets_LOY.sex.txt', sep='')
filename_CNY_diff_PAR1 = paste(output_prefix, '.tmp.facets_LOY.diff_PAR1.txt', sep='')
filename_CNY_diff_PAR2 = paste(output_prefix, '.tmp.facets_LOY.diff_PAR2.txt', sep='')


load(facets_rdata)


df.stat_LOY = df
df.stat_LOY$is_LOX = NA
df.stat_LOY$is_LOY = NA
df.stat_LOY$LOX_median_cncf = NA
df.stat_LOY$LOY_median_cncf = NA
df.stat_LOY$PAR1_additional = NA
df.stat_LOY$PAR2_additional = NA
df.stat_LOY$CNV_normal_X = NA
df.stat_LOY$CNV_tumor_X = NA
df.stat_LOY$CNV_length_X = NA
df.stat_LOY$CNV_length_ratio_X = NA
df.stat_LOY$GAIN_LOSS_X = NA
df.stat_LOY$CNV_normal_Y = NA
df.stat_LOY$CNV_tumor_Y = NA
df.stat_LOY$CNV_length_Y = NA
df.stat_LOY$CNV_length_ratio_Y = NA
df.stat_LOY$GAIN_LOSS_Y = NA
df.stat_LOY$CNV_normal_PAR1 = NA
df.stat_LOY$CNV_tumor_PAR1 = NA
df.stat_LOY$CNV_length_PAR1 = NA
df.stat_LOY$CNV_length_ratio_PAR1 = NA
df.stat_LOY$GAIN_LOSS_PAR1 = NA
df.stat_LOY$diff_PAR1_XY = NA
df.stat_LOY$CNV_normal_PAR2 = NA
df.stat_LOY$CNV_tumor_PAR2 = NA
df.stat_LOY$CNV_length_PAR2 = NA
df.stat_LOY$CNV_length_ratio_PAR2 = NA
df.stat_LOY$GAIN_LOSS_PAR2 = NA
df.stat_LOY$diff_PAR2_XY = NA



df.CNV_XY = cncf %>% dplyr::filter(chrom %in% c(23,24,25,26)) %>% left_join(df_genomesize %>% select(chr_order, chr_num) %>% rename('chrom' = 'chr_order')) %>% mutate(chrom = chr_num)
df.CNV_XY$pair_name = df.stat_LOY$pair_name

tmp.df = df.CNV_XY %>% select(1:6) %>% dplyr::filter(chrom=='X') 
ratio.nhetX = sum(tmp.df$nhet)/sum(tmp.df$num.mark)
tmp.df = df.CNV_XY %>% select(1:6) %>% dplyr::filter(chrom=='Y')
n.site.Y = sum(tmp.df$num.mark)


purity  = df$purity
# if (is.na(purity)) {purity = 1}



if (n.site.Y< 20 ){inferred_sex = 'female'; tmp.normal.Y = 0 }else{inferred_sex = 'male'; tmp.normal.Y = 1}
if (ratio.nhetX<0.01){tmp.normal.X = 1}else{tmp.normal.X = 2 }
print(paste0("inferred sex is :",inferred_sex))


ploidy = df$ploidy
if (is.na(ploidy)) {ploidy = 2}
# if ((round(ploidy)) %%2 == 1){
# ploidy_hapCN_max = (round(ploidy) + 1)/2
# ploidy_hapCN_min  = ploidy_hapCN_max - 1
# }

df.stat_LOY$purity = purity
df.stat_LOY$ploidy = ploidy
# df.stat_LOY$ploidy_hapCN_min = ploidy_hapCN_min
# df.stat_LOY$ploidy_hapCN_max = ploidy_hapCN_max



# female: XX or XO 
# male: XY or XO
# even it's hard to identify, but Y signals sometimes can be observed in males
# male if nhets on X chromosome is less than 1% of num.mark

df.stat_LOY$inferred_sex = inferred_sex


is.odd <- function(x) x %% 2 != 0


tmp.normal.X_ploidy = tmp.normal.X*round(ploidy)
if(is.odd(tmp.normal.X_ploidy))
{ 	CNV_normal_ploidy_min_X = max(1, (tmp.normal.X_ploidy - 1 )/2); 
	CNV_normal_ploidy_max_X = (tmp.normal.X_ploidy + 1)/2
	}else
{
	CNV_normal_ploidy_min_X = (tmp.normal.X_ploidy )/2
	CNV_normal_ploidy_max_X = (tmp.normal.X_ploidy )/2
}



tmp.normal.Y_ploidy = tmp.normal.Y*round(ploidy)
if(is.odd(tmp.normal.Y_ploidy))
{ 	CNV_normal_ploidy_min_Y = max(1, (tmp.normal.Y_ploidy - 1 )/2); 
	CNV_normal_ploidy_max_Y = (tmp.normal.Y_ploidy + 1)/2
	}else
{
	CNV_normal_ploidy_min_Y = (tmp.normal.Y_ploidy )/2
	CNV_normal_ploidy_max_Y = (tmp.normal.Y_ploidy )/2
}

tmp.normal.PAR1 = tmp.normal.X + tmp.normal.Y
tmp.normal.PAR1_ploidy = tmp.normal.PAR1*round(ploidy)
if(is.odd(tmp.normal.PAR1_ploidy))
{ 	CNV_normal_ploidy_min_PAR1 = max(1, (tmp.normal.PAR1_ploidy - 1 )/2); 
	CNV_normal_ploidy_max_PAR1 = (tmp.normal.PAR1_ploidy + 1)/2
	}else
{
	CNV_normal_ploidy_min_PAR1 = (tmp.normal.PAR1_ploidy )/2
	CNV_normal_ploidy_max_PAR1 = (tmp.normal.PAR1_ploidy )/2
}

tmp.normal.PAR2 = tmp.normal.X + tmp.normal.Y

tmp.normal.PAR2_ploidy = tmp.normal.PAR2*round(ploidy)
if(is.odd(tmp.normal.PAR2_ploidy))
{ 	CNV_normal_ploidy_min_PAR2 = max(1, (tmp.normal.PAR2_ploidy - 1 )/2); 
	CNV_normal_ploidy_max_PAR2 = (tmp.normal.PAR2_ploidy + 1)/2
	}else
{
	CNV_normal_ploidy_min_PAR2 = (tmp.normal.PAR2_ploidy )/2
	CNV_normal_ploidy_max_PAR2 = (tmp.normal.PAR2_ploidy )/2
}



df.stat_LOY = df.stat_LOY %>% mutate(CNV_normal_X    = tmp.normal.X)  %>% mutate(CNV_normal_ploidy_min_X = CNV_normal_ploidy_min_X, CNV_normal_ploidy_max_X = CNV_normal_ploidy_max_X)
df.stat_LOY = df.stat_LOY %>% mutate(CNV_normal_Y    = tmp.normal.Y ) %>% mutate(CNV_normal_ploidy_min_Y = CNV_normal_ploidy_min_Y, CNV_normal_ploidy_max_Y = CNV_normal_ploidy_max_Y)
df.stat_LOY = df.stat_LOY %>% mutate(CNV_normal_PAR1 = tmp.normal.PAR1 ) %>% mutate(CNV_normal_ploidy_min_PAR1 = CNV_normal_ploidy_min_PAR1, CNV_normal_ploidy_max_PAR1 = CNV_normal_ploidy_max_PAR1)
df.stat_LOY = df.stat_LOY %>% mutate(CNV_normal_PAR2 = tmp.normal.PAR2 ) %>% mutate(CNV_normal_ploidy_min_PAR2 = CNV_normal_ploidy_min_PAR2, CNV_normal_ploidy_max_PAR2 = CNV_normal_ploidy_max_PAR2)



df.CNV_XY.sample = df.CNV_XY %>% left_join(df.stat_LOY) %>% dplyr::filter(num.mark >=10)






# X
index_sample = 1
df.CNV_XY.sample.X = df.CNV_XY.sample %>% dplyr::filter(chrom=='X')

df.CNV_XY.sample.X.tmp.CNV = df.CNV_XY.sample.X[!(between(df.CNV_XY.sample.X$tcn,CNV_normal_ploidy_min_X, CNV_normal_ploidy_max_X )) , , drop = F] 

df.stat_LOY[index_sample,'CNV_tumor_X'] = df.CNV_XY.sample.X.tmp.CNV %>%  mutate( tmp = paste(chrom, start, end, cnlr.median, tcn, sep = ":"))  %>% .$tmp %>% paste(collapse = '|')
df.stat_LOY[index_sample,'CNV_length_X'] = sum(df.CNV_XY.sample.X.tmp.CNV %>% mutate(length = end - start + 1) %>% .$length)
df.stat_LOY[index_sample,'CNV_length_ratio_X'] = df.stat_LOY[index_sample,'CNV_length_X'] /sum(df.CNV_XY.sample.X %>% mutate(length = end - start + 1) %>% .$length)





n.CNV_XY.sample.X.gain = dim(df.CNV_XY.sample.X %>%  dplyr::filter( CNV_normal_ploidy_max_X <tcn))[1]
n.CNV_XY.sample.X.loss = dim(df.CNV_XY.sample.X %>%  dplyr::filter(CNV_normal_ploidy_min_X >tcn))[1]
if (n.CNV_XY.sample.X.gain   >0 & n.CNV_XY.sample.X.loss == 0){df.stat_LOY[index_sample,'GAIN_LOSS_X'] =  "Gain"}
if (n.CNV_XY.sample.X.gain == 0 & n.CNV_XY.sample.X.loss   >0){df.stat_LOY[index_sample,'GAIN_LOSS_X'] =  "Loss"}
if (n.CNV_XY.sample.X.gain   >0 & n.CNV_XY.sample.X.loss   >0){df.stat_LOY[index_sample,'GAIN_LOSS_X'] =  "Mixed"}
if (n.CNV_XY.sample.X.gain   ==0 & n.CNV_XY.sample.X.loss   ==0){df.stat_LOY[index_sample,'GAIN_LOSS_X'] =  "None"}
if ((df.stat_LOY[index_sample,'GAIN_LOSS_X'] == "Loss") && df.stat_LOY[index_sample,'CNV_length_ratio_X'] >0.99 ){ df.stat_LOY[index_sample,'is_LOX'] = 1; df.stat_LOY[index_sample,'LOX_median_cncf'] = median(rep(df.CNV_XY.sample.X$cnlr.median, df.CNV_XY.sample.X$num.mark))}else( df.stat_LOY[index_sample,'is_LOX'] = 0 )

# Y

df.CNV_XY.sample.Y = df.CNV_XY.sample %>% dplyr::filter(chrom=='Y')

# df.CNV_XY.sample.Y.tmp.CNV = df.CNV_XY.sample.Y %>%  dplyr::filter( !between(tcn,CNV_normal_ploidy_min_Y, CNV_normal_ploidy_max_Y )) 
df.CNV_XY.sample.Y.tmp.CNV = df.CNV_XY.sample.Y[!(between(df.CNV_XY.sample.Y$tcn,CNV_normal_ploidy_min_Y, CNV_normal_ploidy_max_Y )) , , drop = F] 


df.stat_LOY[index_sample,'CNV_tumor_Y'] = df.CNV_XY.sample.Y.tmp.CNV %>%  mutate( tmp = paste(chrom, start, end, cnlr.median, tcn, sep = ":"))  %>% .$tmp %>% paste(collapse = '|')
df.stat_LOY[index_sample,'CNV_length_Y'] = sum(df.CNV_XY.sample.Y.tmp.CNV %>% mutate(length = end - start + 1) %>% .$length)
df.stat_LOY[index_sample,'CNV_length_ratio_Y'] = df.stat_LOY[index_sample,'CNV_length_Y'] /sum(df.CNV_XY.sample.Y %>% mutate(length = end - start + 1) %>% .$length)


n.CNV_XY.sample.Y.gain = dim(df.CNV_XY.sample.Y %>%  dplyr::filter( CNV_normal_ploidy_max_Y <tcn))[1]
n.CNV_XY.sample.Y.loss = dim(df.CNV_XY.sample.Y %>%  dplyr::filter(CNV_normal_ploidy_min_Y >tcn))[1]
if (n.CNV_XY.sample.Y.gain   >0 & n.CNV_XY.sample.Y.loss == 0){df.stat_LOY[index_sample,'GAIN_LOSS_Y'] =  "Gain"}
if (n.CNV_XY.sample.Y.gain == 0 & n.CNV_XY.sample.Y.loss   >0){df.stat_LOY[index_sample,'GAIN_LOSS_Y'] =  "Loss"}
if (n.CNV_XY.sample.Y.gain   >0 & n.CNV_XY.sample.Y.loss   >0){df.stat_LOY[index_sample,'GAIN_LOSS_Y'] =  "Mixed"}
if (n.CNV_XY.sample.Y.gain   ==0 & n.CNV_XY.sample.Y.loss   ==0){df.stat_LOY[index_sample,'GAIN_LOSS_Y'] =  "None"}
if ((df.stat_LOY[index_sample,'GAIN_LOSS_Y'] == "Loss") && df.stat_LOY[index_sample,'CNV_length_ratio_Y'] >0.99 ){ df.stat_LOY[index_sample,'is_LOY'] = 1; df.stat_LOY[index_sample,'LOY_median_cncf'] = median(rep(df.CNV_XY.sample.Y$cnlr.median, df.CNV_XY.sample.Y$num.mark))}else( df.stat_LOY[index_sample,'is_LOY'] = 0 )




# PAR1
df.CNV_XY.sample.PAR1 = df.CNV_XY.sample %>% dplyr::filter(chrom=='XY_PAR1')

# df.CNV_XY.sample.PAR1.tmp.CNV = df.CNV_XY.sample.PAR1 %>%  dplyr::filter( !between(tcn,CNV_normal_ploidy_min_PAR1, CNV_normal_ploidy_max_PAR1 )) 
df.CNV_XY.sample.PAR1.tmp.CNV = df.CNV_XY.sample.PAR1[!between(df.CNV_XY.sample.PAR1$tcn,CNV_normal_ploidy_min_PAR1, CNV_normal_ploidy_max_PAR1 ),,drop = F ]

df.stat_LOY[index_sample,'CNV_tumor_PAR1'] = df.CNV_XY.sample.PAR1.tmp.CNV %>%  mutate( tmp = paste(chrom, start, end, cnlr.median, tcn, sep = ":"))  %>% .$tmp %>% paste(collapse = '|')
df.stat_LOY[index_sample,'CNV_length_PAR1'] = sum(df.CNV_XY.sample.PAR1.tmp.CNV %>% mutate(length = end - start + 1) %>% .$length)
df.stat_LOY[index_sample,'CNV_length_ratio_PAR1'] = df.stat_LOY[index_sample,'CNV_length_PAR1'] /sum(df.CNV_XY.sample.PAR1 %>% mutate(length = end - start + 1) %>% .$length)


n.CNV_XY.sample.PAR1.gain = dim(df.CNV_XY.sample.PAR1 %>%  dplyr::filter( CNV_normal_ploidy_max_PAR1 <tcn))[1]
n.CNV_XY.sample.PAR1.loss = dim(df.CNV_XY.sample.PAR1 %>%  dplyr::filter(CNV_normal_ploidy_min_PAR1 >tcn))[1]
if (n.CNV_XY.sample.PAR1.gain   >0 & n.CNV_XY.sample.PAR1.loss == 0){df.stat_LOY[index_sample,'GAIN_LOSS_PAR1'] =  "Gain"}
if (n.CNV_XY.sample.PAR1.gain == 0 & n.CNV_XY.sample.PAR1.loss   >0){df.stat_LOY[index_sample,'GAIN_LOSS_PAR1'] =  "Loss"}
if (n.CNV_XY.sample.PAR1.gain   >0 & n.CNV_XY.sample.PAR1.loss   >0){df.stat_LOY[index_sample,'GAIN_LOSS_PAR1'] =  "Mixed"}
if (n.CNV_XY.sample.PAR1.gain   ==0 & n.CNV_XY.sample.PAR1.loss   ==0){df.stat_LOY[index_sample,'GAIN_LOSS_PAR1'] =  "None"}









value.CNV_XY.sample.X.near_PAR1 =  df.CNV_XY.sample.X %>% arrange(start)  %>% slice_head( n  = 1) %>% .$cnlr.median
value.CNV_XY.sample.Y.near_PAR1 =  df.CNV_XY.sample.Y %>% arrange(start)  %>% slice_head( n  = 1) %>% .$cnlr.median
value.CNV_XY.sample.PAR1 =  df.CNV_XY.sample.PAR1 %>% arrange(desc(start)) %>% slice_head( n  = 1) %>% .$cnlr.median
if (inferred_sex == 'male') {diff_PAR1_XY =  value.CNV_XY.sample.PAR1 - (log2(2^value.CNV_XY.sample.Y.near_PAR1 + 2^value.CNV_XY.sample.X.near_PAR1) -1)  } else{diff_PAR1_XY = value.CNV_XY.sample.PAR1 - value.CNV_XY.sample.X.near_PAR1 }
if (length(diff_PAR1_XY) > 0 && abs(diff_PAR1_XY) >0.25  ) {df.stat_LOY[index_sample,'PAR1_additional'] = diff_PAR1_XY}
if (length(diff_PAR1_XY) >0 ){df.stat_LOY[index_sample,'diff_PAR1_XY'] = diff_PAR1_XY}






# PAR2
df.CNV_XY.sample.PAR2 = df.CNV_XY.sample %>% dplyr::filter(chrom=='XY_PAR2')

# df.CNV_XY.sample.PAR2.tmp.CNV = df.CNV_XY.sample.PAR2 %>%  dplyr::filter( !between(tcn,CNV_normal_ploidy_min_PAR2, CNV_normal_ploidy_max_PAR2 )) 
df.CNV_XY.sample.PAR2.tmp.CNV = df.CNV_XY.sample.PAR2[!between(df.CNV_XY.sample.PAR2$tcn,CNV_normal_ploidy_min_PAR2, CNV_normal_ploidy_max_PAR2 ),,drop = F ]


df.stat_LOY[index_sample,'CNV_tumor_PAR2'] = df.CNV_XY.sample.PAR2.tmp.CNV %>%  mutate( tmp = paste(chrom, start, end, cnlr.median, tcn, sep = ":"))  %>% .$tmp %>% paste(collapse = '|')
df.stat_LOY[index_sample,'CNV_length_PAR2'] = sum(df.CNV_XY.sample.PAR2.tmp.CNV %>% mutate(length = end - start + 1) %>% .$length)
df.stat_LOY[index_sample,'CNV_length_ratio_PAR2'] = df.stat_LOY[index_sample,'CNV_length_PAR2'] /sum(df.CNV_XY.sample.PAR2 %>% mutate(length = end - start + 1) %>% .$length)


n.CNV_XY.sample.PAR2.gain = dim(df.CNV_XY.sample.PAR2 %>%  dplyr::filter( CNV_normal_ploidy_max_PAR2 <tcn))[1]
n.CNV_XY.sample.PAR2.loss = dim(df.CNV_XY.sample.PAR2 %>%  dplyr::filter(CNV_normal_ploidy_min_PAR2 >tcn))[1]
if (n.CNV_XY.sample.PAR2.gain   >0 & n.CNV_XY.sample.PAR2.loss == 0){df.stat_LOY[index_sample,'GAIN_LOSS_PAR2'] =  "Gain"}
if (n.CNV_XY.sample.PAR2.gain == 0 & n.CNV_XY.sample.PAR2.loss   >0){df.stat_LOY[index_sample,'GAIN_LOSS_PAR2'] =  "Loss"}
if (n.CNV_XY.sample.PAR2.gain   >0 & n.CNV_XY.sample.PAR2.loss   >0){df.stat_LOY[index_sample,'GAIN_LOSS_PAR2'] =  "Mixed"}
if (n.CNV_XY.sample.PAR2.gain   ==0 & n.CNV_XY.sample.PAR2.loss   ==0){df.stat_LOY[index_sample,'GAIN_LOSS_PAR2'] =  "None"}


value.CNV_XY.sample.X.near_PAR2 =  df.CNV_XY.sample.X %>% arrange(desc(start))  %>% slice_head( n  = 1) %>% .$cnlr.median
value.CNV_XY.sample.Y.near_PAR2 =  df.CNV_XY.sample.Y %>% arrange(desc(start))  %>% slice_head( n  = 1) %>% .$cnlr.median
value.CNV_XY.sample.PAR2 =  df.CNV_XY.sample.PAR2 %>% arrange(start) %>% slice_head( n  = 1) %>% .$cnlr.median
if (inferred_sex == 'male') {diff_PAR2_XY =  value.CNV_XY.sample.PAR2 - (log2(2^value.CNV_XY.sample.Y.near_PAR2 + 2^value.CNV_XY.sample.X.near_PAR2) -1)  } else{diff_PAR2_XY = value.CNV_XY.sample.PAR2 - value.CNV_XY.sample.X.near_PAR2 }

if (length(diff_PAR2_XY) > 0 && abs(diff_PAR2_XY) >0.375) {df.stat_LOY[index_sample,'PAR2_additional'] = diff_PAR2_XY}
if (length(diff_PAR2_XY) >0 ){df.stat_LOY[index_sample,'diff_PAR2_XY'] = diff_PAR2_XY}



## write the files
df.stat_LOY %>% write_tsv(filename_CNY_stat)
write(df.stat_LOY$is_LOX,   filename_CNY_is_LOX)
write(df.stat_LOY$is_LOY,   filename_CNY_is_LOY)
write(df.stat_LOY$inferred_sex,   filename_CNY_inferred_sex)
write(df.stat_LOY$diff_PAR1,   filename_CNY_diff_PAR1)
write(df.stat_LOY$diff_PAR2,   filename_CNY_diff_PAR2)


