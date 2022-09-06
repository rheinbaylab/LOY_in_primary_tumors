library(ggplot2)
library(devtools)
library(stringi)
library(pctGCdata)
library(tidyverse)
# load_all('facetschrY4')
library('facetschrY4')



args = commandArgs(trailingOnly=TRUE)
pair_name = args[1]
pileup = args[2]
min_normal_depth = as.numeric(args[3])
cval = as.numeric(args[4])
maxiter = as.numeric(args[5])
seed_initial = as.numeric(args[6])
seed_iterations = as.numeric(args[7])
print(seed_iterations)

genome_build = "hg19_Y"
load(system.file("ext_genomesize", "genomesize.rda", package = "facetschrY4"))
load(system.file("ext_genomesize", "hg19_Y.rda", package = "facetschrY4"))
df_genomesize = list.genomesize[[genome_build]]
df_genomesize = df_genomesize %>% mutate(if_auto = if_else(chr_order <=22, 1, 0))

# pair_name='test-TP-NB'
# pileup='test-TP-NB.pileup.add_mid_exon.gz
# min_normal_depth = 20
# cval = 150
# maxiter = 1
# seed_initial = 10
# seed_iterations = 1



if (seed_iterations <= 0) {
    seed_iterations = 1
}

run_facets = function(seed, pileup, min_normal_depth, cval, genome_build, maxiter) {
    set.seed(seed)
    genome_build = "hg19_Y"
    load(system.file("ext_genomesize", "genomesize.rda", package = "facetschrY4"))
    load(system.file("ext_genomesize", "hg19_Y.rda", package = "facetschrY4"))
    df_genomesize = list.genomesize[[genome_build]]
df_genomesize = df_genomesize %>% mutate(if_auto = if_else(chr_order <=22, 1, 0))
    # rcmat = readSnpMatrix(pileup)
    # rcmat = f_get_split_chr(rcmat, df.par, vector.chr_order)
    # xx = preProcSample(rcmat, ndepth = min_normal_depth, cval = cval, gbuild = "udef", ugcpct  = hg19_Y)
    # oo = procSample(xx, cval = cval)
    # fit = emcncf(oo, maxiter = maxiter)

  rcmat = readSnpMatrix(pileup)
  rcmat = f_get_split_chr(rcmat, df.par, vector.chr_order)
  rcmat = rcmat %>% left_join(df_genomesize %>% dplyr::select(chr_num,chr_order), by = c('Chromosome' = 'chr_num')) %>% mutate(Chromosome = chr_order)
  rcmat = rcmat[,1:6]
  list.out_pre = preProcSample(rcmat, ndepth = min_normal_depth, cval = cval, gbuild = "udef", ugcpct  = hg19_Y)
  oo = procSample(list.out_pre, cval = cval)
  fit = emcncf(oo, maxiter = maxiter)
  out = oo

    return(list(seed_oo=oo, seed_fit=fit))
}

plot_facets_iterations = function(pair_name, seeds_dataframe, median_purity, median_ploidy) {
    title = paste('Stability of FACETS outputs across ', seed_iterations, ' seeds,\n', pair_name, sep='')
    p = ggplot(seeds_dataframe, aes(x=ploidy, y=purity)) + xlim(0, 8) + ylim(0, 1) +
        geom_hline(yintercept=median_purity, linetype='dashed', size=1, alpha=0.2) +
        geom_vline(xintercept=median_ploidy, linetype='dashed', size=1, alpha=0.2) +
        geom_jitter(size = 6, color='#E69F00', alpha=0.5) +
        ggtitle(title) +
        theme(plot.title = element_text(size=16, hjust=0.5)) +
        theme(axis.title.x = element_text(size=16)) +
        theme(axis.title.y = element_text(size=16))

    return(p)
}



seeds = seed_initial:(seed_initial+seed_iterations-1)
seeds_dataframe = data.frame()
for (seed in seeds) {
    print(seed)
    seed_list = run_facets(seed, pileup, min_normal_depth, cval, genome_build, maxiter)
    seed_oo = seed_list$seed_oo
    seed_fit = seed_list$seed_fit

    dip_log_r = seed_oo$dipLogR
    purity = as.numeric(signif(seed_fit$purity, 3))
    ploidy = as.numeric(signif(seed_fit$ploidy, 3))

    seed_dataframe = data.frame(seed, purity, ploidy, dip_log_r)
    names(seed_dataframe) <- c("seed", "purity", "ploidy", "dip_log_r")
    seeds_dataframe = rbind(seeds_dataframe, seed_dataframe)
}


seeds_dataframe_filename = paste(pair_name, '.facets_iterations.txt', sep='')
write.table(seeds_dataframe, seeds_dataframe_filename, sep='\t', quote=FALSE, row.names=FALSE)

idx_na_purity = is.na(seeds_dataframe$purity)
n_iteratons_na_purity = sum(as.numeric(idx_na_purity))
median_purity = median(seeds_dataframe[!idx_na_purity, ]$purity, na.rm=TRUE)
median_ploidy = median(seeds_dataframe[!idx_na_purity, ]$ploidy, na.rm=TRUE)

plot_iterations_filename = paste(pair_name, '.facets_iterations.pdf', sep='')
pdf(plot_iterations_filename, height=10, width=7.5)
plot_facets_iterations(pair_name, seeds_dataframe, median_purity, median_ploidy)
dev.off()

if (!is.na(median_purity)) {
  delta_median_purity = abs(seeds_dataframe$purity - median_purity)
  min_delta_purity = min(delta_median_purity, na.rm=TRUE)
  median_purity_seed = as.numeric(seeds_dataframe[delta_median_purity %in% min_delta_purity,]$seed)
  if (length(median_purity_seed) > 1) {median_purity_seed = median_purity_seed[1]}
  used_seed = as.numeric(median_purity_seed)
} else {
  used_seed = seed_initial
}



facets_list = run_facets(used_seed, pileup, min_normal_depth, cval, genome_build, maxiter)
oo = facets_list$seed_oo
fit = facets_list$seed_fit

purity = as.numeric(signif(fit$purity, 3))
ploidy = as.numeric(signif(fit$ploidy, 3))
log_likelihood = fit$loglik
cncf = fit$cncf
cncf = cncf %>% left_join(oo$out %>% dplyr::select(seg, cf, tcn,lcn))
dip_log_r = oo$dipLogR
flags = oo$flags
emflags = fit$emflags

number_segments = NROW(cncf)
number_segments_NA_LCN = sum(is.na(cncf$lcn.em))

genome_segments_filename_auto = paste(pair_name, '.genome_segments.facets.auto.pdf', sep='')
genome_segments_filename_all = paste(pair_name, '.genome_segments.facets.allchr.pdf', sep='')
diagnostic_plot_filename = paste(pair_name, '.diagnostic_plot.pdf', sep='')
cncf_dataframe_filename = paste(pair_name, '.facets_cncf.tsv', sep='')
summary_dataframe_filename = paste(pair_name, '.facets_summary.tsv', sep='')
flags_filename = paste(pair_name, '.tmp.facets_flags.txt', sep='')
emflags_filename = paste(pair_name, '.tmp.facets_emflags.txt', sep='')


if(!is.na(purity)) {
    # pdf.options(reset = TRUE, onefile = FALSE)
    # pdf(genome_segments_filename_auto, height=10, width=7.5)
    # plotSample(x=oo, emfit=fit)
    # dev.off()

    pdf.options(reset = TRUE, onefile = FALSE)
    pdf(genome_segments_filename_auto, height=10, width=7.5)
    oo_new = oo
    fit_new = fit
    oo_new$jointseg = oo_new$jointseg %>% dplyr::filter(chrom %in% c(1:23))
    oo_new$out = oo_new$out %>% dplyr::filter(chrom %in% c(1:23))
    oo_new$chromlevels = c(1:23)
    fit_new$cncf = fit_new$cncf %>% dplyr::filter(chrom %in% c(1:23))
    plotSample(x=oo_new, emfit=fit_new)
    dev.off()


    p = plotSample_chrsize(x = oo, emfit = fit, df_genomesize = df_genomesize)
    ggsave(p, filename =genome_segments_filename_all, height = 10, width = 7.5, dpi = 400)


    pdf(diagnostic_plot_filename, height=10, width=7.5)
    logRlogORspider(oo$out, oo$dipLogR)
    dev.off()
} else {
    write.table(NA, genome_segments_filename_auto)
    write.table(NA, genome_segments_filename_all)
    write.table(NA, diagnostic_plot_filename)
}



df = data.frame(pair_name=pair_name, purity=purity, ploidy=ploidy, digLogR=dip_log_r, log_likelihood = ifelse(is.null(log_likelihood), "", log_likelihood), flags = paste(flags,collapse = '|'), emflags = paste(emflags,collapse = '|'), number_segments=number_segments, number_segments_NA_LCN=number_segments_NA_LCN)


df_flags = data.frame(flags=flags)
df_emflags = data.frame(emflags=emflags)


write.table(cncf, cncf_dataframe_filename, sep='\t', quote=FALSE, row.names=FALSE)
write.table(df, summary_dataframe_filename, sep='\t', quote=FALSE, row.names=FALSE)
write.table(df_flags, flags_filename, sep='\t', quote=FALSE, row.names=FALSE)
write.table(df_emflags, emflags_filename, sep='\t', quote=FALSE, row.names=FALSE)

write(purity, paste0(pair_name,'.tmp.purity.txt'))
write(ploidy, paste0(pair_name, '.tmp.ploidy.txt'))
write(log_likelihood, paste0(pair_name, '.tmp.log_likelihood.txt'))
write(dip_log_r, paste0(pair_name, '.tmp.dip_log_r.txt'))
write(used_seed, paste0(pair_name, '.tmp.seed_used.txt'))
write(n_iteratons_na_purity, paste0(pair_name, '.tmp.number_iterations_with_na_purity.txt'))
write(number_segments, paste0(pair_name, '.tmp.number_segments.txt'))
write(number_segments_NA_LCN, paste0(pair_name, '.tmp.number_segments_NA_LCN.txt'))


# save RData object needed for add_ccf_to_maf_config method, part of Phylogic preprocessing
formatSegmentOutput = function(out, sampID) {
    seg = list()
    seg$ID = rep(sampID, nrow(out$out))
    seg$chrom = out$out$chrom
    seg$loc.start = rep(NA, length(seg$ID))
    seg$loc.end = seg$loc.start
    seg$num.mark = out$out$num.mark
    seg$seg.mean = out$out$cnlr.median
    for(i in 1:nrow(out$out)) {
        lims = range(out$jointseg$maploc[(out$jointseg$chrom == out$out$chrom[i] & out$jointseg$seg == out$out$seg[i])], na.rm=T)
        seg$loc.start[i] = lims[1]
        seg$loc.end[i] = lims[2]
    }
    as.data.frame(seg)
}

out = oo
out$IGV <- formatSegmentOutput(out, pair_name)
save(fit, out,df, cncf,  file=paste(pair_name, '.facets.RData', sep=''))

