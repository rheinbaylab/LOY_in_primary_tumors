plotSample_ggplot <-  function(x, emfit=NULL, clustered=FALSE, plot.type=c("em","naive","both","none"), sname=NULL){



    jseg <- x$jointseg
    chrbdry <- which(diff(jseg$chrom) != 0)

        out <- emfit$cncf
        # add the naive tcn, lcn and cf to out
        out$tcn <- x$out$tcn
        out$lcn <- x$out$lcn
        out$cf <- x$out$cf

    if (clustered) {
        cnlr.median <- out$cnlr.median.clust
        mafR <- out$mafR.clust
        mafR[is.na(mafR)] <- out$mafR[is.na(mafR)]
    } else {
        cnlr.median <- out$cnlr.median
        mafR <- out$mafR
    }
    mafR <- abs(mafR)
    # chromosome colors
    chrcol <- 1+rep(out$chrom-2*floor(out$chrom/2), out$num.mark)
    nn <- cumsum(table(jseg$chrom[is.finite(jseg$cnlr)]))
    segbdry <- cumsum(c(0,out$num.mark))
    segstart <- segbdry[-length(segbdry)]
    segend <- segbdry[-1]



color_deep = c("lightblue","grey","azure4","slateblue")

theme_custom = function(){
  theme_bw() +
  theme(panel.spacing = unit(0, "lines")) +
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  theme(strip.background = element_blank()) +
  theme(plot.background = element_blank()) +
  theme(strip.text = element_blank()) + 
  theme(panel.grid = element_blank()) + 
  theme(axis.text.x  = element_blank()) + 
  theme(axis.ticks.x = element_blank())
}

jseg$chrcol = rep((1:x$nX+1)%%2 +1 , times = as.vector(table(jseg$chrom))) 
jseg$chrcol = factor(jseg$chrcol, levels = c(1,2))
# jseg %>% mutate(order = 1:nrow(jseg))  %>% ggplot(aes(x = order, y = cnlr,color = chrcol) )+ geom_point()  + scale_colour_manual(values = color_deep) 

# jseg = jseg %>% dplyr::filter(chrom %in% c(23,24))
jseg = jseg %>% mutate(order = 1:nrow(jseg))

# p1 = jseg %>% mutate(order = 1:nrow(jseg))  %>% ggplot(aes(x = order, y = cnlr,color = chrcol) )+ geom_point()  + scale_colour_manual(values = color_deep) + facet_grid(cols = vars(chrom),scale = 'free_x',space = 'free_x') + theme_custom()


df_segment = data.frame(segstart, cnlr.median, segend)

p1 = jseg  %>% ggplot(aes(x = order, y = cnlr,color = chrcol) )+ geom_point()  + scale_colour_manual(values = color_deep) + theme_custom() + 
xlab('') + ylab('log-ratio') + theme(legend.position = 'none')+ 
geom_hline(yintercept = median(jseg$cnlr, na.rm=TRUE), colour="green2") + 
geom_vline(xintercept = chrbdry,lwd = 0.25) + geom_hline(yintercept = x$dipLogR, col = "magenta4")+ 
geom_segment(data = df_segment, aes(x = segstart, y = cnlr.median, xend = segend, yend = cnlr.median), colour =  2,lwd = 2)

    cfpalette <- c(colorRampPalette(c("white", "steelblue"))(10),"bisque2")

        cfcol <- cfpalette[round(10*out$cf+0.501)]
df_segment1 = data.frame(segstart,segend)
df_segment1$lcn.em = out$lcn.em
df_segment1$tcn.em = out$tcn.em
df_segment1$lcn = out$lcn
df_segment1$tcn = out$tcn
df_segment1$plus_mafR = sqrt(mafR)
df_segment1$minus_mafR = -1 * sqrt(mafR)
df_segment1 = df_segment1 %>% mutate(xmin = segstart, xmax = segend, ymin = 0, ymax = 1 )
df_segment1$cfcol = cfcol

p2 = jseg  %>% ggplot(aes(x = order, y = valor,color = chrcol) )+ geom_point()  + scale_colour_manual(values = color_deep) + theme_custom() + 
xlab('') + ylab('log-odds-ratio') + theme(legend.position = 'none')+ 
geom_hline(yintercept = median(jseg$cnlr, na.rm=TRUE), colour="green2") + 
geom_vline(xintercept = chrbdry,lwd = 0.25) + geom_hline(yintercept = x$dipLogR, col = "magenta4")+ 
geom_segment(data = df_segment1, aes(x = segstart, y = plus_mafR, xend = segend, yend = plus_mafR), colour =  2,lwd = 2) + 
geom_segment(data = df_segment1, aes(x = segstart, y = minus_mafR, xend = segend, yend = minus_mafR), colour =  2,lwd = 2) 






p3 = ggplot() + geom_segment(data = df_segment1, aes(x = segstart, y = lcn, xend = segend, yend = lcn), colour =  2,lwd = 2) + 
geom_segment(data = df_segment1, aes(x = segstart, y = tcn, xend = segend, yend = tcn), colour =  1,lwd = 2) + 
geom_vline(xintercept = chrbdry, lwd = 0.25) + ylab('copy number (em)') + theme_custom() + xlab('')




theme_custom_p4 = function(){
  theme_bw() +
  theme(panel.spacing = unit(0, "lines")) +
  theme(strip.background = element_blank()) +
  theme(plot.background = element_blank()) +
  theme(strip.text = element_blank()) + 
  theme(panel.grid = element_blank()) 
}

    chromlevels <- x$chromlevels

    chromlevels_label = chromlevels %>% as.data.frame() %>% rename('chr_order' = '.') %>% left_join(df_genomesize) %>% .$chr_num

p4 = ggplot(df_segment1,aes(xmin = xmin, xmax = xmax , ymin = ymin, ymax = ymax,fill = cfcol)) + geom_rect() + theme_custom_p4()  + theme(legend.position = 'none') + ylab('cf_em') + xlab('Chromosome') + scale_fill_manual(values = cfcol)  + coord_cartesian(ylim = c(0,1), clip = 'off') + scale_y_continuous(expand = c(0,0)) + scale_x_continuous(labels=chromlevels_label,breaks = (nn+c(0,nn[-length(nn)]))/2) + theme(axis.title.y = element_text(angle = 0, hjust = 1, vjust = 0.5), axis.text.y = element_blank(),axis.ticks.y = element_blank())


# grid.arrange(p1, p2,  p3, p4,ncol=1, nrow = 4,  heights = c(9,9,6,3))
plot_grid(p1, p2, p3, p4,ncol=1, nrow = 4, rel_heights = c(9,9,6,3),align = 'v', axis = 'lr')
}










plotSample_chrsize <- function(x, emfit=NULL, clustered=FALSE, plot.type=c("em","naive","both","none"), sname=NULL, df_genomesize = df_genomesize)
{

    # x = oo
    # emfit = fit
    # clustered = FALSE
    # plot.type = 'em'
    # sname = NULL
    # df_genomesize = df_genomesize

      chromlevels <- x$chromlevels
    chromlevels_label = chromlevels %>% as.data.frame() %>% rename('chr_order' = '.') %>% left_join(df_genomesize) %>% .$chr_num
    jseg <- x$jointseg
    chrbdry =c(0, df_genomesize$chr_cumsum)
    df_genomesize$x_chrlength  = c(0,df_genomesize$chr_cumsum[-c(nrow(df_genomesize))])

    out <- emfit$cncf
    # add the naive tcn, lcn and cf to out
    out$tcn <- x$out$tcn
    out$lcn <- x$out$lcn
    out$cf <- x$out$cf

    cnlr.median <- out$cnlr.median
    mafR <- out$mafR
    mafR <- abs(mafR)
    # chromosome colors
    chrcol <- 1+rep(out$chrom-2*floor(out$chrom/2), out$num.mark)
    nn <- cumsum(table(jseg$chrom[is.finite(jseg$cnlr)]))

    out = out %>% left_join(df_genomesize, by = c('chrom' = 'chr_order')) %>% mutate(start_x = start + x_chrlength, end_x = end+ x_chrlength)
    segstart = out$start_x
    segend = out$end_x




      color_deep = c("lightblue","grey","azure4","slateblue")

      theme_custom = function(){
        theme_bw() +
        theme(panel.spacing = unit(0, "lines")) +
        theme(plot.margin = unit(c(0,0,0,0), "lines")) +
        theme(strip.background = element_blank()) +
        theme(plot.background = element_blank()) +
        theme(strip.text = element_blank()) + 
        theme(panel.grid = element_blank()) + 
        theme(axis.text.x  = element_blank()) + 
        theme(axis.ticks.x = element_blank())
      }
      n_chr = length(as.vector(table(jseg$chrom)))
      # jseg$chrcol = rep((1:x$nX+1)%%2 +1 , times = as.vector(table(jseg$chrom))) 
      jseg$chrcol = rep((1:n_chr+1)%%2 +1 , times = as.vector(table(jseg$chrom))) 
      jseg$chrcol = factor(jseg$chrcol, levels = c(1,2))

      jseg = jseg %>% mutate(order = 1:nrow(jseg))


      df_segment = data.frame(segstart, cnlr.median, segend)
      jseg = jseg %>% left_join(df_genomesize, by = c('chrom' = 'chr_order')) %>% mutate(maploc_x = maploc + x_chrlength)

      p1 = jseg %>%  ggplot(aes(x = maploc_x, y = cnlr,color = chrcol) )+ geom_point()  + scale_colour_manual(values = color_deep) + theme_custom() + 
      xlab('') + ylab('log-ratio') + theme(legend.position = 'none')+ 
      geom_hline(yintercept = median(jseg$cnlr, na.rm=TRUE), colour="green2") + 
      geom_vline(xintercept = chrbdry,lwd = 0.25) + geom_hline(yintercept = x$dipLogR, col = "magenta4")+ 
      geom_segment(data = df_segment, aes(x = segstart, y = cnlr.median, xend = segend, yend = cnlr.median), colour =  2,lwd = 2)

      cfpalette <- c(colorRampPalette(c("white", "steelblue"))(10),"bisque2")

      cfcol <- cfpalette[round(10*out$cf+0.501)]
      df_segment1 = data.frame(segstart,segend)
      df_segment1$lcn.em = out$lcn.em
      df_segment1$tcn.em = out$tcn.em
      df_segment1$lcn = out$lcn
      df_segment1$tcn = out$tcn
      df_segment1$plus_mafR = sqrt(mafR)
      df_segment1$minus_mafR = -1 * sqrt(mafR)
      df_segment1 = df_segment1 %>% mutate(xmin = segstart, xmax = segend, ymin = 0, ymax = 1 )
      df_segment1$cfcol = cfcol

      p2 = jseg  %>% ggplot(aes(x = maploc_x, y = valor,color = chrcol) )+ geom_point()  + scale_colour_manual(values = color_deep) + theme_custom() + 
      xlab('') + ylab('log-odds-ratio') + theme(legend.position = 'none')+ 
      geom_hline(yintercept = median(jseg$cnlr, na.rm=TRUE), colour="green2") + 
      geom_vline(xintercept = chrbdry,lwd = 0.25) + geom_hline(yintercept = x$dipLogR, col = "magenta4")+ 
      geom_segment(data = df_segment1, aes(x = segstart, y = plus_mafR, xend = segend, yend = plus_mafR), colour =  2,lwd = 2) + 
      geom_segment(data = df_segment1, aes(x = segstart, y = minus_mafR, xend = segend, yend = minus_mafR), colour =  2,lwd = 2) 



      p3 = ggplot() + geom_segment(data = df_segment1, aes(x = segstart, y = lcn, xend = segend, yend = lcn), colour =  2,lwd = 2) + 
      geom_segment(data = df_segment1, aes(x = segstart, y = tcn, xend = segend, yend = tcn), colour =  1,lwd = 2) + 
      geom_vline(xintercept = chrbdry, lwd = 0.25) + ylab('copy number (em)') + theme_custom() + xlab('')




      theme_custom_p4 = function(){
        theme_bw() +
        theme(panel.spacing = unit(0, "lines")) +
        theme(strip.background = element_blank()) +
        theme(plot.background = element_blank()) +
        theme(strip.text = element_blank()) + 
        theme(panel.grid = element_blank()) 
      }



      p4 = ggplot(df_segment1,aes(xmin = xmin, xmax = xmax , ymin = ymin, ymax = ymax,fill = cfcol)) + geom_rect() + theme_custom_p4()  + theme(legend.position = 'none') + ylab('cf_em') + xlab('Chromosome') + scale_fill_manual(values = cfcol)  + coord_cartesian(ylim = c(0,1), clip = 'off') + scale_y_continuous(expand = c(0,0)) + scale_x_continuous(labels=df_genomesize$chr_num,breaks = df_genomesize$mid_point) + theme(axis.title.y = element_text(angle = 0, hjust = 1, vjust = 0.5), axis.text.y = element_blank(),axis.ticks.y = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + xlab('')


      # grid.arrange(p1, p2,  p3, p4,ncol=1, nrow = 4,  heights = c(9,9,6,3))
      plot_grid(p1, p2, p3, p4,ncol=1, nrow = 4, rel_heights = c(9,9,6,3),align = 'v', axis = 'lr') + theme(panel.background = element_rect(fill = 'white',color = 'white'))
}












plotSample_chrsize_PAR <- function(x, emfit=NULL, clustered=FALSE, plot.type=c("em","naive","both","none"), sname=NULL, df_genomesize = df_genomesize, chr_selected = vector.selected_chr)
{

    # x = oo
    # emfit = fit
    # clustered = FALSE
    # plot.type = 'em'
    # sname = NULL
    # df_genomesize = df_genomesize
    # chr_selected = vector.selected_chr

    chromlevels <- x$chromlevels
    chromlevels_label = chromlevels %>% as.data.frame() %>% rename('chr_order' = '.') %>% left_join(df_genomesize) %>% .$chr_num
    jseg <- x$jointseg
    chrbdry =c(0, df_genomesize$chr_cumsum)
    df_genomesize$x_chrlength  = c(0,df_genomesize$chr_cumsum[-c(nrow(df_genomesize))])

    out <- emfit$cncf
    # add the naive tcn, lcn and cf to out
    out$tcn <- x$out$tcn
    out$lcn <- x$out$lcn
    out$cf <- x$out$cf

    cnlr.median <- out$cnlr.median
    mafR <- out$mafR
    mafR <- abs(mafR)
    # chromosome colors
    chrcol <- 1+rep(out$chrom-2*floor(out$chrom/2), out$num.mark)
    nn <- cumsum(table(jseg$chrom[is.finite(jseg$cnlr)]))

    out = out %>% left_join(df_genomesize, by = c('chrom' = 'chr_order')) %>% mutate(start_x = start + x_chrlength, end_x = end+ x_chrlength)
    segstart = out$start_x
    segend = out$end_x




    color_deep = c("lightblue","grey","azure4","slateblue")

    theme_custom = function(){
      theme_bw() +
      theme(panel.spacing = unit(0, "lines")) +
      theme(plot.margin = unit(c(0,0,0,0), "lines")) +
      theme(strip.background = element_blank()) +
      theme(plot.background = element_blank()) +
      theme(strip.text = element_blank()) + 
      theme(panel.grid = element_blank()) + 
      theme(axis.text.x  = element_blank()) + 
      theme(axis.ticks.x = element_blank())
    }

    jseg$chrcol = rep((1:x$nX+1)%%2 +1 , times = as.vector(table(jseg$chrom))) 
    jseg$chrcol = factor(jseg$chrcol, levels = c(1,2))

    jseg = jseg %>% mutate(order = 1:nrow(jseg))


    df_segment = data.frame(segstart, cnlr.median, segend)


df_genomesize_selected = df_genomesize %>% dplyr::filter(chr_num %in% chr_selected)
df_genomesize_selected = df_genomesize_selected %>% mutate(length = ifelse(grepl('PAR',chr),5000000,length)) %>% mutate(chr_cumsum = cumsum(as.numeric(length)))
    df_genomesize_selected$start_chrlength  = c(0,df_genomesize_selected$chr_cumsum[-c(nrow(df_genomesize_selected))])
    df_genomesize_selected = df_genomesize_selected %>% mutate(mid_point  = chr_cumsum - length/2)
    jseg = jseg %>% left_join(df_genomesize_selected, by = c('chrom' = 'chr_order')) %>% mutate(maploc_x = maploc + x_chrlength) %>% dplyr::filter(chr_num %in% chr_selected) 


    p1 = jseg %>% dplyr::filter(chr_num %in% chr_selected) %>%  ggplot(aes(x = maploc_x, y = cnlr,color = chrcol) )+ geom_point()  + scale_colour_manual(values = color_deep) + theme_custom() + 
    xlab('') + ylab('log-ratio') + theme(legend.position = 'none')+ 
    geom_hline(yintercept = median(jseg$cnlr, na.rm=TRUE), colour="green2") + 
    geom_vline(xintercept = chrbdry,lwd = 0.25) + geom_hline(yintercept = x$dipLogR, col = "magenta4")+ 
    geom_segment(data = df_segment, aes(x = segstart, y = cnlr.median, xend = segend, yend = cnlr.median), colour =  2,lwd = 2) %>% 



    cfpalette <- c(colorRampPalette(c("white", "steelblue"))(10),"bisque2")

    cfcol <- cfpalette[round(10*out$cf+0.501)]
    df_segment1 = data.frame(segstart,segend)
    df_segment1$lcn.em = out$lcn.em
    df_segment1$tcn.em = out$tcn.em
    df_segment1$lcn = out$lcn
    df_segment1$tcn = out$tcn
    df_segment1$plus_mafR = sqrt(mafR)
    df_segment1$minus_mafR = -1 * sqrt(mafR)
    df_segment1 = df_segment1 %>% mutate(xmin = segstart, xmax = segend, ymin = 0, ymax = 1 )
    df_segment1$cfcol = cfcol

    p2 = jseg  %>% ggplot(aes(x = maploc_x, y = valor,color = chrcol) )+ geom_point()  + scale_colour_manual(values = color_deep) + theme_custom() + 
    xlab('') + ylab('log-odds-ratio') + theme(legend.position = 'none')+ 
    geom_hline(yintercept = median(jseg$cnlr, na.rm=TRUE), colour="green2") + 
    geom_vline(xintercept = chrbdry,lwd = 0.25) + geom_hline(yintercept = x$dipLogR, col = "magenta4")+ 
    geom_segment(data = df_segment1, aes(x = segstart, y = plus_mafR, xend = segend, yend = plus_mafR), colour =  2,lwd = 2) + 
    geom_segment(data = df_segment1, aes(x = segstart, y = minus_mafR, xend = segend, yend = minus_mafR), colour =  2,lwd = 2) 



    p3 = ggplot() + geom_segment(data = df_segment1, aes(x = segstart, y = lcn, xend = segend, yend = lcn), colour =  2,lwd = 2) + 
    geom_segment(data = df_segment1, aes(x = segstart, y = tcn, xend = segend, yend = tcn), colour =  1,lwd = 2) + 
    geom_vline(xintercept = chrbdry, lwd = 0.25) + ylab('copy number (em)') + theme_custom() + xlab('')




    theme_custom_p4 = function(){
      theme_bw() +
      theme(panel.spacing = unit(0, "lines")) +
      theme(strip.background = element_blank()) +
      theme(plot.background = element_blank()) +
      theme(strip.text = element_blank()) + 
      theme(panel.grid = element_blank()) 
    }



    p4 = ggplot(df_segment1,aes(xmin = xmin, xmax = xmax , ymin = ymin, ymax = ymax,fill = cfcol)) + geom_rect() + theme_custom_p4()  + theme(legend.position = 'none') + ylab('cf_em') + xlab('Chromosome') + scale_fill_manual(values = cfcol)  + coord_cartesian(ylim = c(0,1), clip = 'off') + scale_y_continuous(expand = c(0,0)) + scale_x_continuous(labels=chromlevels_label,breaks = df_genomesize$mid_point) + theme(axis.title.y = element_text(angle = 0, hjust = 1, vjust = 0.5), axis.text.y = element_blank(),axis.ticks.y = element_blank())


    # grid.arrange(p1, p2,  p3, p4,ncol=1, nrow = 4,  heights = c(9,9,6,3))
    plot_grid(p1, p2, p3, p4,ncol=1, nrow = 4, rel_heights = c(9,9,6,3),align = 'v', axis = 'lr')
}
