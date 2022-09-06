plotSample_chrsize_single_chromosome <- function(list.out, input_chr = 23, emfit=NULL, clustered=FALSE, plot.type="both", sname=NULL, df_genomesize = df_genomesize, y_axis = TRUE)
{

    # list.out = oo
    # emfit = fit
    # clustered = FALSE
    # plot.type = 'em'
    # sname = NULL
    # df_genomesize = df_genomesize
    # input_chr = 23
    # y_axis = FALSE

    # set for chromosome
    chromlevels <- list.out$chromlevels
    chromlevels_label = chromlevels %>% as.data.frame() %>% rename('chr_order' = '.') %>% left_join(df_genomesize) %>% .$chr_num
    df_genomesize_single = df_genomesize %>% dplyr::filter(chr_order == !!input_chr)
    # chromosome names
    chr_name = df_genomesize_single %>% dplyr::filter(chr_order == !!input_chr ) %>% .$chr_num

    # set for snp position
    df.snp <- list.out$jointseg %>% dplyr::filter(chrom == !!input_chr)
    df.snp = df.snp  %>% mutate(chrcol = chrom %%2 + 1) 
    df.snp$chrcol = factor(df.snp$chrcol, levels = c(1,2))
    if (nrow(df.snp) >0) {df.snp = df.snp %>% arrange(maploc) %>% mutate(order = 1:nrow(df.snp))}else{df.snp <-  df.snp %>% add_column(order = NA)}
    x_chrlength = 0
    df.snp = df.snp %>% left_join(df_genomesize_single, by = c('chrom' = 'chr_order')) %>% mutate(maploc_x = maploc + !!x_chrlength)
    # set the value for extreme data
    df.snp <- df.snp %>% mutate(cnlr = if_else(abs(cnlr)>10, sign(cnlr)*10, cnlr))
    df.snp <- df.snp %>% mutate(valor = if_else(abs(valor)>10, sign(valor)*10, valor))
    df.snp <- df.snp %>% mutate(lorvar = if_else(abs(lorvar)>10, sign(lorvar)*10, lorvar))

    # add the em/naive tcn, lcn and cf to out
    df.seg <- emfit$cncf
    df.seg <- df.seg %>% dplyr::left_join(list.out$out %>% dplyr::select(seg,cf, tcn, lcn))
    df.seg$abs_mafR  = abs(df.seg$mafR)
    # set the value for extreme data
    df.seg <- df.seg %>% mutate(cnlr.median = if_else(abs(cnlr.median)>10, sign(cnlr.median)*10, cnlr.median))
    df.seg <- df.seg %>% mutate(cf.em = if_else(abs(cf.em)>10, sign(cf.em)*10, cf.em))
    df.seg <- df.seg %>% mutate(tcn.em = if_else(abs(tcn.em)>10, sign(tcn.em)*10, tcn.em))
    df.seg <- df.seg %>% mutate(lcn.em = if_else(abs(lcn.em)>10, sign(lcn.em)*10, lcn.em))
    df.seg <- df.seg %>% mutate(cf = if_else(abs(cf)>10, sign(cf)*10, cf))
    df.seg <- df.seg %>% mutate(tcn = if_else(abs(tcn)>10, sign(tcn)*10, tcn))
    df.seg <- df.seg %>% mutate(lcn = if_else(abs(lcn)>10, sign(lcn)*10, lcn))
    df.seg <- df.seg %>% mutate(abs_mafR = if_else(abs(abs_mafR)>10, sign(abs_mafR)*10, abs_mafR))
    # select for the chromosome data
    df.seg = df.seg %>% dplyr::filter(chrom==!!input_chr)
    df.seg = df.seg %>% left_join(df_genomesize_single, by = c('chrom' = 'chr_order')) %>% mutate(start_x = start + !!x_chrlength, end_x = end+ !!x_chrlength)

    # set for the cooridnate limit 
    tmp.max_cnlr = min(max(abs(list.out$jointseg %>% dplyr::filter(!is.na(cnlr))  %>% .$cnlr)), 10)
    tmp.max_valor = min(max(abs(list.out$jointseg %>% dplyr::filter(!is.na(valor)) %>% .$valor)),10)
    tmp.max_lcn_tcn = min(max(list.out$out$tcn,list.out$out$lcn,emfit$cncf$tcn.em, emfit$cncf$lcn.em,na.rm = T),10)
    total.median.cnlr = list.out$jointseg$cnlr
    total.median.dipLogR = list.out$dipLogR

    color_deep = c("lightblue","grey","azure4","slateblue")

    theme_custom = function(y_axis){
      theme_y_axis =  theme_bw() +
      theme(panel.spacing = unit(0, "lines")) +
      theme(plot.margin = unit(c(0,0,0,0), "lines")) +
      theme(strip.background = element_blank()) +
      theme(plot.background = element_blank()) +
      theme(strip.text = element_blank()) + 
      theme(panel.grid = element_blank()) + 
      theme(axis.text.x  = element_blank()) + 
      theme(axis.ticks.x = element_blank()) 

      if(y_axis){
        return(theme_y_axis)
      }else
      {return(theme_y_axis +  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) )}
    }


    # abline(v=chrbdry, lwd=0.25)
    # abline(h=median(jseg$cnlr, na.rm=TRUE), col="green2")
    # abline(h = x$dipLogR, col = "magenta4")

  p1 = df.snp %>%  ggplot(aes(x = maploc_x, y = cnlr,color = chrcol) )+ geom_point()  + scale_colour_manual(values = color_deep) + theme_custom(y_axis) + 
  xlab('') + ylab('log-ratio') + theme(legend.position = 'none') + 
  geom_hline(yintercept = median(total.median.cnlr, na.rm=TRUE), colour="green2") + geom_hline(yintercept = total.median.dipLogR, colour = 'black', linetype = 'dashed')  + geom_hline(yintercept = c(1,-1), colour = 'grey', linetype = 'dashed') + 
  geom_segment(data = df.seg, aes(x = start_x, y = cnlr.median, xend = end_x, yend = cnlr.median), colour =  2,lwd = 2) + ylim(-tmp.max_cnlr, tmp.max_cnlr)


  p2 = df.snp  %>% ggplot(aes(x = maploc_x, y = valor,color = chrcol) )+ geom_point()  + scale_colour_manual(values = color_deep) + theme_custom(y_axis) + 
  xlab('') + ylab('log-odds-ratio') + theme(legend.position = 'none')+ 
  geom_hline(yintercept = median(total.median.cnlr, na.rm=TRUE), colour="green2") + geom_hline(yintercept = total.median.dipLogR, colour = 'black', linetype = 'dashed')  +geom_hline(yintercept = c(1,-1), colour = 'grey', linetype = 'dashed') + 
  geom_segment(data = df.seg , aes(x = start_x, y = sqrt(abs_mafR), xend = end_x, yend = sqrt(abs_mafR)), colour =  2,lwd = 2) + 
  geom_segment(data = df.seg , aes(x = start_x, y = -sqrt(abs_mafR), xend = end_x, yend = -sqrt(abs_mafR)), colour =  2,lwd = 2) + 
  ylim(-tmp.max_valor, tmp.max_valor)

  p3 = ggplot() + geom_segment(data = df.seg, aes(x = start_x, y = lcn.em, xend = end_x, yend = lcn.em), colour =  2,lwd = 2) + 
  geom_segment(data = df.seg, aes(x = start_x, y = tcn.em, xend = end_x, yend = tcn.em), colour =  1,lwd = 2) + 
  ylab('copy number (em)') + theme_custom(y_axis) + xlab('') + xlab(chr_name) + ylim(-0.5, tmp.max_lcn_tcn)

  p4 = ggplot() + geom_segment(data = df.seg, aes(x = start_x, y = lcn, xend = end_x, yend = lcn), colour =  2,lwd = 2) + 
  geom_segment(data = df.seg, aes(x = start_x, y = tcn, xend = end_x, yend = tcn), colour =  1,lwd = 2) + 
  ylab('copy number') + theme_custom(y_axis) + xlab('') + xlab(chr_name) + ylim(-0.5, tmp.max_lcn_tcn)


  p5 = plot_grid(p1, p2, p3, p4, ncol=1, nrow = 4, rel_heights = c(9,9,9,9),align = 'v', axis = 'lr')
  p5
}




library(cowplot)
library(tidyverse)
library(ggplot2)
library(pctGCdata)
library(facetschrY4)
library(devtools)

# load_all('facetschrY4')
library('facetschrY4')


genome_build = "hg19_Y"
load(system.file("ext_genomesize", "genomesize.rda", package = "facetschrY4"))
load(system.file("ext_genomesize", "hg19_Y.rda", package = "facetschrY4"))
df_genomesize = list.genomesize[[genome_build]]


args = commandArgs(trailingOnly=TRUE)

facets_out = args[1] # can be output from facets or from pileup
output_prefix = args[2]
sample_title = args[3]



fig_XY_output_pdf = paste0(output_prefix,'.facets.CNY.chrXY.pdf', sep = '')
fig_XY_output_png = paste0(output_prefix,'.facets.CNY.chrXY.png', sep = '')
fig_all_output_png = paste0(output_prefix,'.facets.CNY.all_chr.png', sep = '')
output_facets_rdata = paste0(output_prefix,'.RData', sep = '')



if(endsWith(facets_out,'RData') ){
  load(facets_out)
  list.out = out
}else{
  pileup = facets_out
  rcmat = readSnpMatrix(pileup)
  rcmat = f_get_split_chr(rcmat, df.par, vector.chr_order)
  rcmat = rcmat %>% left_join(df_genomesize %>% dplyr::select(chr_num,chr_order), by = c('Chromosome' = 'chr_num')) %>% mutate(Chromosome = chr_order)
  rcmat = rcmat[,1:6]
  # rcmat_new = rcmat %>% dplyr::filter(Chromosome %in% c(4,25))
  list.out_pre = preProcSample(rcmat, ndepth = 20, cval = 250, gbuild = "udef", ugcpct  = hg19_Y)
  list.out = procSample(list.out_pre, cval = 250)
  fit = emcncf(list.out, maxiter = 5)
  # save the data and plot for the original facets plot figure
  out = list.out
  save(fit, out,list.out , file = output_facets_rdata)

  # # plot for original figure
  # genome_segments_filename = paste0(output_prefix,'.raw_facts.pdf',sep = '')
  # pdf(genome_segments_filename, height=10, width=7.5)
  # plotSample(x=list.out, emfit=fit)
  # dev.off()

}


# plot for sex chromosomes
pX = plotSample_chrsize_single_chromosome(list.out = list.out, emfit = fit, df_genomesize = df_genomesize ,input_chr = 23)
pY = plotSample_chrsize_single_chromosome(list.out = list.out, emfit = fit, df_genomesize = df_genomesize ,input_chr = 24, y_axis = FALSE) 
p_par1 = plotSample_chrsize_single_chromosome(list.out = list.out, emfit = fit, df_genomesize = df_genomesize ,input_chr = 25, y_axis = FALSE)
p_par2 = plotSample_chrsize_single_chromosome(list.out = list.out, emfit = fit, df_genomesize = df_genomesize ,input_chr = 26, y_axis = FALSE)

p4 = plot_grid(pX, pY, p_par1,p_par2 ,nrow = 1, rel_widths = c(8,8,4,4)) + theme(panel.background = element_rect(fill = 'white',color = 'white'))


  title <- ggdraw() + 
    draw_label(
      sample_title,
      x = 0,
      hjust = 0
    ) +
    theme(
      plot.margin = margin(0, 0, 0, 0)
      # plot.margin = margin(0, 0, 0, 7)
    )

  p5 = plot_grid(
    title, p4,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 1)
  ) + theme(panel.background = element_rect(fill = 'white',color = 'white'))




ggsave(p5, filename = fig_XY_output_pdf, dpi = 300,width = 10, height = 7)
ggsave(p5, filename = fig_XY_output_png, dpi = 300,width = 10, height = 7)



# plot for other chromosome
p1 = plotSample_chrsize_single_chromosome(list.out = list.out , emfit = fit, df_genomesize = df_genomesize, input_chr = 1)
p2 = plotSample_chrsize_single_chromosome(list.out = list.out , emfit = fit, df_genomesize = df_genomesize, input_chr = 2, y_axis = FALSE )
p3 = plotSample_chrsize_single_chromosome(list.out = list.out , emfit = fit, df_genomesize = df_genomesize, input_chr = 3, y_axis = FALSE )
p4 = plotSample_chrsize_single_chromosome(list.out = list.out , emfit = fit, df_genomesize = df_genomesize, input_chr = 4, y_axis = FALSE )
p5 = plotSample_chrsize_single_chromosome(list.out = list.out , emfit = fit, df_genomesize = df_genomesize, input_chr = 5, y_axis = FALSE )
p6 = plotSample_chrsize_single_chromosome(list.out = list.out , emfit = fit, df_genomesize = df_genomesize, input_chr = 6, y_axis = FALSE )
p7 = plotSample_chrsize_single_chromosome(list.out = list.out , emfit = fit, df_genomesize = df_genomesize, input_chr = 7, y_axis = FALSE )
p8 = plotSample_chrsize_single_chromosome(list.out = list.out , emfit = fit, df_genomesize = df_genomesize, input_chr = 8, y_axis = FALSE )
p9 = plotSample_chrsize_single_chromosome(list.out = list.out , emfit = fit, df_genomesize = df_genomesize, input_chr = 9, y_axis = FALSE )
p10 = plotSample_chrsize_single_chromosome(list.out = list.out , emfit = fit, df_genomesize = df_genomesize, input_chr = 10, y_axis = FALSE )
p11 = plotSample_chrsize_single_chromosome(list.out = list.out , emfit = fit, df_genomesize = df_genomesize, input_chr = 11, y_axis = FALSE )
p12 = plotSample_chrsize_single_chromosome(list.out = list.out , emfit = fit, df_genomesize = df_genomesize, input_chr = 12, y_axis = FALSE )
p13 = plotSample_chrsize_single_chromosome(list.out = list.out , emfit = fit, df_genomesize = df_genomesize, input_chr = 13, y_axis = FALSE )
p14 = plotSample_chrsize_single_chromosome(list.out = list.out , emfit = fit, df_genomesize = df_genomesize, input_chr = 14, y_axis = FALSE )
p15 = plotSample_chrsize_single_chromosome(list.out = list.out , emfit = fit, df_genomesize = df_genomesize, input_chr = 15, y_axis = FALSE )
p16 = plotSample_chrsize_single_chromosome(list.out = list.out , emfit = fit, df_genomesize = df_genomesize, input_chr = 16, y_axis = FALSE )
p17 = plotSample_chrsize_single_chromosome(list.out = list.out , emfit = fit, df_genomesize = df_genomesize, input_chr = 17, y_axis = FALSE )
p18 = plotSample_chrsize_single_chromosome(list.out = list.out , emfit = fit, df_genomesize = df_genomesize, input_chr = 18, y_axis = FALSE )
p19 = plotSample_chrsize_single_chromosome(list.out = list.out , emfit = fit, df_genomesize = df_genomesize, input_chr = 19, y_axis = FALSE )
p20 = plotSample_chrsize_single_chromosome(list.out = list.out , emfit = fit, df_genomesize = df_genomesize, input_chr = 20, y_axis = FALSE )
p21 = plotSample_chrsize_single_chromosome(list.out = list.out , emfit = fit, df_genomesize = df_genomesize, input_chr = 21, y_axis = FALSE )
p22 = plotSample_chrsize_single_chromosome(list.out = list.out , emfit = fit, df_genomesize = df_genomesize, input_chr = 22, y_axis = FALSE )
pX = plotSample_chrsize_single_chromosome(list.out = list.out, emfit = fit, df_genomesize = df_genomesize ,input_chr = 23, y_axis = FALSE)
pY = plotSample_chrsize_single_chromosome(list.out = list.out, emfit = fit, df_genomesize = df_genomesize ,input_chr = 24, y_axis = FALSE) 
p_par1 = plotSample_chrsize_single_chromosome(list.out = list.out, emfit = fit, df_genomesize = df_genomesize ,input_chr = 25, y_axis = FALSE)
p_par2 = plotSample_chrsize_single_chromosome(list.out = list.out, emfit = fit, df_genomesize = df_genomesize ,input_chr = 26, y_axis = FALSE)


p_all = plot_grid(p1 , p2 , p3 , p4 , p5 , p6 , p7 , p8 , p9 , p10 , p11 , p12 , p13 , p14 , p15 , p16 , p17 , p18 , p19 , p20 , p21 , p22,pX,pY, p_par1, p_par2 ,nrow = 1, rel_widths = c(2,rep(1, 25))) + theme(panel.background = element_rect(fill = 'white',color = 'white'))


  title <- ggdraw() + 
    draw_label(
      sample_title,
      x = 0,
      hjust = 0
    ) +
    theme(
      plot.margin = margin(0, 0, 0, 0)
      # plot.margin = margin(0, 0, 0, 7)
    )

  p_all_with_title = plot_grid(
    title, p_all,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 1)
  ) + theme(panel.background = element_rect(fill = 'white',color = 'white'))



ggsave(p_all_with_title, filename = fig_all_output_png, dpi = 300,width = 27, height = 7)
