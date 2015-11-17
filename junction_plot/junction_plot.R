require(gridExtra)
require(ggplot2)
require(reshape2)
require(Hmisc)

make_differential_splicing_plot=function(y, x, len=500, length_transform=function(g) log(g+1), main_title=NA, snp_pos=NA) {

  introns=colnames(y)
  intron_meta=do.call(rbind,strsplit(introns,":"))
  colnames(intron_meta)=c("chr","start","end","clu")
  intron_meta=as.data.frame(intron_meta,stringsAsFactors = F)
  intron_meta$start=as.numeric(intron_meta$start)
  intron_meta$end=as.numeric(intron_meta$end)
  
  #jxn_names=paste0(intron_meta$chr,":",intron_meta$start,"-",intron_meta$end)
  #jxn_counts=colSums(y)
  
  #names(jxn_counts)=NULL
  #require(rPython)
  #python.load("lib/sashimi_plot_utils.py")
  #python.call("plot_junctions_forR", jxn_names, jxn_counts)

  new_theme_empty <- theme_bw(base_size = 16)
  new_theme_empty$line <- element_blank()
  new_theme_empty$rect <- element_blank()
  new_theme_empty$strip.text <- element_blank()
  new_theme_empty$axis.text <- element_blank()
  
                                        #new_theme_empty$axis.title <- element_blank()

  groups=sort(unique(x))
  
  max_log=.5*ceil(2*log10( 1+max( unlist( foreach (tis=groups) %do% { intron_meta$counts=colMeans(y[ tis==x,,drop=F]) } ) ) ))

  breaks=if (max_log <= 2.5) seq(0,max_log,by=0.5) else seq(0,ceil(max_log),by=1)
  limits=c(0.0,max_log)

  intron_meta$id=as.factor(1:nrow(intron_meta))
  temp=intron_meta[,c("id","start","end")]
  m=melt(temp, id.vars = "id")

  s=unique(m$value)
  if (!is.na(snp_pos)) s=c(s,snp_pos)
  s=sort(s)
  d=s[2:length(s)]-s[1:length(s)-1]
  trans_d=length_transform(d) # e.g. log(d+1), sqrt(d), d ^ .33
  coords=c(0,cumsum(trans_d)) 
  names(coords)=s

  snp_coord=coords[as.character(snp_pos)]

  total_length=sum(trans_d) # ==max(coords)

  first_plot=T
  
  min_height=0
  max_height=0
  
  plots=foreach (tis=groups) %do% {
    print(tis)
    intron_meta$counts=colMeans(y[ tis==x,,drop=F])
    group_sample_size=sum(tis==x)
    print(intron_meta$counts)
    
    allEdges=do.call(rbind,foreach (i=1:nrow(intron_meta)) %do% {
      start=coords[ as.character(intron_meta$start[i]) ]
      end=coords[ as.character(intron_meta$end[i]) ]
      l=end-start
      #edge = data.frame(bezier(x=c(start, (start + end)/2, end), y=c(0, (1+sqrt(l)) * ( (i %% 2)*2-1 ), 0),evaluation = len))
      h=(1+sqrt(l)) * ( (i %% 2)*2-1 )
      min_height=min(min_height,h)
      max_height=max(max_height,h)
      edge = data.frame(bezier(x=c(start, start-0.1*total_length, (start + end)/2, end+0.1*total_length, end), y=c(0, h*2/3, h, h*2/3, 0),evaluation = len))
      edge$Sequence <- log10(1+intron_meta$counts[i]) * sin( seq(0,pi,length.out=len) ) # For size and colour weighting in plot
      edge$log10counts=log10(intron_meta$counts[i])
      edge$Group <- i
      edge
    })

    if (is.na(main_title) | !first_plot) new_theme_empty$plot.title <- element_blank()
    first_plot=F
    
    g=ggplot(allEdges) + geom_path(aes(x = x, y = y, group = Group, colour=log10counts, size = Sequence, alpha=.9)) + scale_size(breaks=breaks, labels=format(10^breaks,digits=0), limits=limits, range = c(.3, 10), guide = guide_legend(title = "mean count")) + scale_alpha(guide="none",range = c(0.1, 1)) + new_theme_empty + scale_color_gradient(breaks=breaks, limits=limits, labels=format(10^breaks,digits=0),low="blue",high="red", guide = guide_legend(title = "mean count")) + ylab(paste0(tis," (n=",group_sample_size,")")) + xlab("") + xlim(-.05*total_length,1.05*total_length) + geom_hline(yintercept=0,alpha=.3) #  + scale_color_discrete(guide="none")
    if (!is.na(snp_coord)) {
        df=data.frame(x=snp_coord,xend=snp_coord,y=0,yend=max_height*1.1)
        print(df)
        g=g + geom_segment(data=df,aes(x=x,y=y,xend=xend,yend=yend)) #+geom_vline(xintercept=snp_coord)
    }
    g
  }
  
  df=data.frame(x=coords, xend=total_length*(s-min(s))/(max(s)-min(s)), y=0, yend=min_height)
  plots[[length(plots)]]=plots[[length(plots)]] + geom_segment(data=df, aes(x=x,y=y,xend=xend,yend=yend),alpha=.3)
  
  exons_chr=exons[exons==intron_meta$chr[1],]
  exons_here=exons_chr[ ( min(s) <= exons_chr$start & exons_chr$start <= max(s) ) | ( min(s) <= exons_chr$end & exons_chr$end <= max(s) ), ]

  df=data.frame(x=total_length*(exons_here$start-min(s))/(max(s)-min(s)), xend=total_length*(exons_here$end-min(s))/(max(s)-min(s)), y=min_height, yend=min_height)
  #invert_mapping=function(pos) if (pos %in% s) coords[as.character(pos)] else {
  #  w=which( pos < s[2:length(s)] & pos > s[1:(length(s)-1)] )
  #  coords[w] + (coords[w+1]-coords[w])*(pos - s[w])/(s[w+1]-s[w])
  #}
  cat(dim(exons_here))
  plots[[length(plots)]]=plots[[length(plots)]] + geom_segment(data=df, aes(x=x,y=y,xend=xend,yend=yend), alpha=.3, size=5)  + geom_hline(yintercept=min_height,alpha=.3)
  
  if (!is.na(main_title)) plots[[1]] = plots[[1]] + ggtitle(main_title)
  
  do.call( grid.arrange, c(plots, list(ncol=1)))
}



