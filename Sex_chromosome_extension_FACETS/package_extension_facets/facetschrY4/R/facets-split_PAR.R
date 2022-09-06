f_get_split_chr = function(rcmat, df.par, vector.chr_order = vector.chr_order)
{
rcmat.new = rcmat %>% mutate(Chromosome =
 ifelse(Chromosome == 'X' & between(Position, df.par$start[1],df.par$end[1]) ,'XY_PAR1',Chromosome)) %>% mutate(Chromosome = ifelse(Chromosome =='X' & between(Position, df.par$start[2],df.par$end[2]), 'XY_PAR2', Chromosome)) %>% 
	dplyr::filter(!(Chromosome == 'Y' &  (between(Position, df.par$start[3],df.par$end[3]) | between(
		Position, df.par$start[4],df.par$end[4])))) 

rcmat.new = rcmat.new %>% dplyr::filter(Chromosome %in% vector.chr_order )
rcmat.new$Chromosome = factor(rcmat.new$Chromosome, levels = vector.chr_order)
rcmat.new = rcmat.new %>% arrange(Chromosome)
return(rcmat.new)

}
