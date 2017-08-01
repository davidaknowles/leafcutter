
#' Make cluster level plot
#'
#' @import ggplot2
#' @import gridExtra
#' @export 
make_cluster_plot <- function(
  cluster_to_plot, 
  main_title = NA, 
  exons_table = NULL, 
  meta = NULL, 
  cluster_ids = NULL, 
  counts = NULL, 
  introns = NULL,
  snp_pos=NA){
  
  if( is.null(cluster_to_plot)){
    print("no cluster selected!")
    grid.arrange(tableGrob("hello there"))
  }
  
  meta$group=as.factor(meta$group)
  group_names=levels(meta$group)
  
  stopifnot(cluster_to_plot %in% cluster_ids)
  # create variables for later
  y <- t(counts[ cluster_ids==cluster_to_plot, ])
  #x <- numeric(nrow(y))+1
  x <- meta$group
  length_transform <- function(g){ log(g+1) }
  introns_to_plot <- introns[ introns$clusterID == cluster_to_plot, ] 
  summary_func <- colSums
  legend_title <- "Mean counts"
  alphabet <- c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z")
  
  
  junction_colour <- "red"
  cryptic_colour <- "pink"
  
  print(introns_to_plot)
  
  # convert colnames(y) into meta data
  intron_meta=get_intron_meta(colnames(y))
  
  # if introns_to_plot is present, add the verdict
  if(!is.null(introns_to_plot)){
    intron_meta$verdict <- introns_to_plot$verdict[match(paste(intron_meta$start,intron_meta$end), paste(introns_to_plot$start, introns_to_plot$end) ) ]
    intron_meta$dPSI <- introns_to_plot$deltapsi[match(paste(intron_meta$start,intron_meta$end), paste(introns_to_plot$start, introns_to_plot$end) ) ]
  # give alphabetical rank based on dPSI
    ranks <- alphabet[1:nrow(intron_meta)]
    absPSI <- intron_meta$dPSI[ order( abs(intron_meta$dPSI), decreasing=TRUE) ]
    #ranks <- cbind(ranks,absPSI)
    #intron_meta$rank <- ranks[ order(abs(intron_meta$dPSI), decreasing=TRUE) ]
    intron_meta$rank <- ranks[match(intron_meta$dPSI, absPSI)]
  }
  print(intron_meta)
  #return(intron_meta)
  new_theme_empty <- theme_bw(base_size = 15 )
  new_theme_empty$panel.background = element_rect(fill="white", colour = "white")
  new_theme_empty$line <- element_blank()
  new_theme_empty$rect <- element_blank()
  new_theme_empty$strip.text <- element_blank()
  new_theme_empty$axis.text <- element_blank()
  
  groups=sort(unique(x))
  
  max_log=.5*ceiling(2*log10( 1+max( unlist( foreach (tis=groups) %do% { intron_meta$counts=summary_func(y[ tis==x,,drop=F]) } ) ) ))
  
  breaks=if (max_log <= 2.5) seq(0,max_log,by=0.5) else seq(0,ceiling(max_log),by=1)
  limits=c(0.0,max_log)
  
  intron_meta$id=as.factor(1:nrow(intron_meta)) # number each junction
  temp=intron_meta[,c("id","start","end")]
  m=reshape2::melt(temp, id.vars = "id") # melt to get list of all coordinates
  
  s=unique(m$value) # get unique start and end values
  if (!is.na(snp_pos)) s=c(s,snp_pos)
  s=sort(s) 
  d=s[2:length(s)]-s[1:length(s)-1] # get the difference between pairs of coordinates
  trans_d <- length_transform(d) # e.g. log(d+1), sqrt(d), d ^ .33 # apply a trasnformation function - doesn't work on negative numbers!
  coords <- c(0,cumsum(trans_d))
  names(coords)=s
  
  snp_coord=coords[as.character(snp_pos)]
  
  total_length=sum(trans_d) # ==max(coords)
  my_xlim=c(-.05*total_length,1.05*total_length)
  first_plot=T
 
  #print(intron_meta)
  ############## 
  # PLOT SETTINGS
  ###############
  
  min_height=0
  max_height=0
  curv <- 0.1
  min_exon_length <- 0.5
  maxratio=0
  minratio=1.0
  yFactor = 0.65   # originally set as 0.65
  yConstant = -0.25 # originally set as 0.5
  labelTextSize=3.5 # orignally set as 5
  curveMax = 10
  curveExponent = 2
  yOffset = 0
  # horizontal white line to clean up the edges of the curves
  centreLineWidth = 3
  
  
  # for gene name colouring (black by default, other colours for multiple genes)
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  #summary_func=function(a) apply(a,2,mean)
  #summary_func=function(a) apply(sweep(a,1,rowSums(a),"/"),2,mean)
  # sweep is dividing each entry in each row by the sum of all entries in that row and then apply is finding the mean value of each column
  summary_func=function(a) apply( sweep(a,1,rowSums(a),"/"),2, function(g) mean(g, na.rm=T) ) 

  plots <- foreach (tis=groups) %do% {
    # print(y[tiss=x,,drop=F])
    intron_meta$counts=summary_func(y[ tis==x,,drop=F])
    #print(length(intron_meta$counts))
    maxratio=max(c(max(intron_meta$counts/sum(intron_meta$counts)), maxratio)) # find max and min values of summarised counts
    minratio=min(c(min(intron_meta$counts/sum(intron_meta$counts)), minratio))
  }
  #print(c(maxratio, minratio))
  last_group=groups[length(groups)]
  
  plots <- list()
  for( fancyVar in 1:length(groups) ){

    intron_meta$counts=summary_func(y[ groups[fancyVar]==x,,drop=F])
    intron_meta$prop=intron_meta$counts#/sum(intron_meta$counts)  # this has been changed

    group_sample_size=sum(groups[fancyVar]==x)
    #print(intron_meta$counts)

    # (l^yFactor/2) + 0.5


    #print(coords)

    # for each junction
    allEdges=do.call(rbind,foreach (i=1:nrow(intron_meta)) %do% {
      #allEdges=do.call(rbind,foreach (i=1:2) %do% {
      if (i%%2==1) return(NULL)  # only care about the even numbered junctions?
      #if (intron_meta$counts[i]==0) return(NULL)
      start=coords[ as.character(intron_meta$start[i]) ]
      end=coords[ as.character(intron_meta$end[i]) ]
      l=end-start  # length of junction
      #edge = data.frame(bezier(x=c(start, (start + end)/2, end), y=c(0, (1+sqrt(l)) * ( (i %% 2)*2-1 ), 0),evaluation = len))
      # h=(1+sqrt(l)) * ( (i %% 2)*2-1 )
      # min_height=min(min_height,h)
      # max_height=max(max_height,h)

      edge = data.frame(start, end)
      edge$startv <- intron_meta$start[i]
      edge$endv <- intron_meta$end[i]
      edge$start <- start  # why do this again? You already set them
      edge$end <- end
      edge$log10counts=intron_meta$counts[i]+1

      # label each junction
      edge$label=paste0(format(intron_meta$prop[i],digits=2, scientific=FALSE),"^", intron_meta$rank[i])
      #edge$label=format(intron_meta$counts[i], scientific=F)
      #edge$label=paste(format(intron_meta$counts[i], scientific=F),'\n(',format(intron_meta$prop[i],digits=2),')',sep='')
      edge$clu<- intron_meta$clu[i]
      #edge$Curve <- j
      edge$Group <- i
      edge$xtext <-start+l/2
      edge$ytext <- -( ( l^(yFactor) / 2 ) + yConstant)  # magic formula here
      #edge$SIZE <- intron_meta$prop[i]
      #edge$SIZE <- intron_meta$prop[i]+1 # make SIZE equal to the normalised count + 1
      #edge$SIZE <- as.factor(.bincode(intron_meta$prop[i]+0.1, breaks=seq(0,100,1)/100, include.lowest=TRUE))
      edge$verdict <- ifelse( intron_meta$verdict[i] == "annotated", yes = "annotated", no ="cryptic")
      edge
    })

    allEdgesP=do.call(rbind,foreach (i=1:nrow(intron_meta)) %do% {
      if (i%%2==0) return(NULL)  # just the odd numbered junctions
      #if (intron_meta$counts[i]==0) return(NULL)
      start=coords[ as.character(intron_meta$start[i]) ]
      end=coords[ as.character(intron_meta$end[i]) ]
      l=end-start
      #edge = data.frame(bezier(x=c(start, (start + end)/2, end), y=c(0, (1+sqrt(l)) * ( (i %% 2)*2-1 ), 0),evaluation = len))
      # h=(1+sqrt(l)) * ( (i %% 2)*2-1 )
      # min_height=min(min_height,h)
      # max_height=max(max_height,h)

      edge = data.frame(start, end)
      edge$startv <- intron_meta$start[i]
      edge$endv <- intron_meta$end[i]
      edge$start <- start
      edge$end <- end
      edge$log10counts=intron_meta$counts[i]+1
      edge$label=paste0(format(intron_meta$prop[i],digits=2, scientific=FALSE),"^",intron_meta$rank[i])
      #edge$label=format(intron_meta$prop[i],digits=2, scientific=T)
      #edge$label=paste(format(intron_meta$counts[i], scientific=F),'\n(',format(intron_meta$prop[i],digits=2),')',sep='')
      edge$clu<- intron_meta$clu[i]
      #edge$Curve <- j
      edge$Group <- i
      edge$xtext <-start+l/2
      edge$ytext <- l^(yFactor)/2+yConstant
      #edge$SIZE <- intron_meta$prop[i]
      #edge$SIZE <- as.factor(.bincode(intron_meta$prop[i], breaks=seq(0,100,1)/100, include.lowest=TRUE))
      edge$SIZE <- intron_meta$prop[i]+1
      edge$verdict <- ifelse( intron_meta$verdict[i] == "annotated", yes = "annotated", no ="cryptic")
      edge
    })
    
    #print(allEdges)
    #print(allEdgesP)
    #print(max(c(max(allEdges$SIZE), max(allEdgesP$SIZE))))
    if (is.na(main_title) | !first_plot) new_theme_empty$plot.title <- element_blank()
    first_plot=F

    # YLIMN=8 # set as 8
    # YLIMP=-9 # set as -9

    MAXcounts <- max(c(max(with(allEdgesP,1)),max(with(allEdges,1))))
    #YLIMN <- max(with(allEdgesP,end-start))
    #YLIMP <- max(with(allEdges,end-start))
    #print(is.numeric(allEdgesP$SIZE))
    #print(allEdges)
    #print(allEdgesP)

    YLIMP <- max( allEdgesP$ytext) + 0.25 * max( allEdgesP$ytext)
    YLIMN <- min( allEdges$ytext) + 0.25 * min( allEdges$ytext)

    # print(YLIMN)
    # print(YLIMP)
    #
    #
    # print(allEdgesP)
    # print(allEdges)
    # if value is 0 then don't plot curve
    #allEdgesP <- allEdgesP[ allEdgesP$label != 0, ]
    #allEdges <- allEdges[ allEdges$label != 0, ]
    # junction_colour if annotated? a different colour if not
    g <- ggplot() +
      geom_curve(data=allEdgesP, aes(x = start, xend = xtext, y = yOffset, yend = ytext, group = Group, colour = verdict, size = curveMax * (log10counts - 1)^curveExponent ),
                 angle=90, curvature=-curv,lineend="round") +
      geom_curve(data=allEdgesP, aes(x = xtext, xend = end, y = ytext, yend = yOffset, group = Group, colour = verdict, size = curveMax * (log10counts - 1)^curveExponent ),
                 angle=90, curvature=-curv,lineend="round") +
      geom_curve(data=allEdges, aes(x = start, xend = xtext, y = -yOffset, yend = ytext, group = Group, colour = verdict, size = curveMax * (log10counts - 1)^curveExponent ),
                 angle=90,curvature=curv,lineend="round") +
      geom_curve(data=allEdges, aes(x = xtext, xend = end, y = ytext, yend = -yOffset, group = Group, colour = verdict, size = curveMax * (log10counts - 1)^curveExponent ),
                 angle=90,curvature=curv,lineend="round") +
      
      

      #geom_curve(data=allEdges,aes(x = xtext+floor(log10(log10counts))/2, xend = end, y = ytext, yend = 0, group = Group,  color=SIZE), size =2,
      #		angle=40,curvature=-curv,lineend="round",ncp=1000) +

      # label the start and end coordinates above the line
      #geom_text(data=allEdgesP,aes(x=start,y=5,label=startv,angle=90, size =3))+
      #geom_text(data=allEdges,aes(x=end,y=-5,label=endv, size=3, angle = 90))+
      # label the start and end coordinates below the line
      #geom_text(data=allEdges,aes(x=start,y=-5,label=startv,angle=90, size =3))+
      # geom_text(data=allEdgesP,aes(x=end,y=5,label=endv, size=3, angle = 90))+

    # put some vertical lines to mark the start and end points
    #geom_linerange(data=allEdgesP,aes(x=start,ymin = 0, ymax=2.5),alpha = 0.2)+
    #geom_linerange(data=allEdgesP,aes(x=end,ymin = 0, ymax=2.5), alpha = 0.2)+
    #geom_linerange(data=allEdges,aes(x=start,ymin = -2.5, ymax=0),alpha = 0.2)+fm
    #geom_linerange(data=allEdges,aes(x=end,ymin = -2.5, ymax=0),alpha = 0.2)+
    new_theme_empty +
      # make the y axis label the group
      ylab(paste0(groups[fancyVar]," (n=",group_sample_size,")")) +
      xlab("") +
      xlim(my_xlim) +
      # try titling instead - why doesn't this work?
      ggtitle(paste0(groups[fancyVar]," (n=",group_sample_size,")" ) ) +

      # horizontal line - smooth out the ends of the curves
      geom_hline(yintercept=0, size = centreLineWidth, colour = "white") +
      geom_hline(yintercept=0,alpha=.9, size=1) +

      # label the junctions
      geom_label(data=allEdgesP,aes(x=xtext,y=0.95*ytext,label=label), size = labelTextSize, label.size = 0, parse=TRUE, fill = "white",colour = "black", label.r = unit(0.3,"lines"), label.padding = unit(0.3,"lines") ) +
      geom_label(data=allEdges,aes(x=xtext,y=0.95*ytext,label=label), size= labelTextSize, label.size = 0, parse=TRUE, fill = "white", colour = "black", label.r = unit(0.3,"lines"), label.padding = unit(0.3,"lines") ) +
      #
      ylim(YLIMN,YLIMP) +
      scale_size_continuous(limits=c(0,10),guide='none')
    # is this used for anything? color is currently set to clu which doesn't change for each junction
    #scale_color_discrete( c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928'),guide='none') #+ scale_color_gradient2(low='white',mid='lightblue',high='purple4',guide="none")  #  + scale_color_discrete(guide="none")

    #g=ggplot(allEdges) + geom_path(aes(x = x, y = y, group = Group, colour=log10counts, size = Sequence, alpha=.9)) + scale_size(breaks=breaks, labels=format(10^breaks,digits=0), limits=limits, range = c(.3, 10), guide = guide_legend(title = legend_title)) + scale_alpha(guide="none",range = c(0.1, 1)) + new_theme_empty + scale_color_gradient(breaks=breaks, limits=limits, labels=format(10^breaks,digits=0),low="orange2",high="darkred", guide = guide_legend(title = legend_title)) + ylab(paste0(tis," (n=",group_sample_size,")")) + xlab("") + xlim(my_xlim) + geom_hline(yintercept=0,alpha=.3) #  + scale_color_discrete(guide="none")
    if (!is.na(snp_coord)) {
      df=data.frame(x=snp_coord,xend=snp_coord,y=0,yend=max_height*1.1)
      #print(df)
      g <- g + geom_segment(data=df,aes(x=x,y=y,xend=xend,yend=yend)) #+geom_vline(xintercept=snp_coord)
    }
    plots[[fancyVar]] <- g  # return the plot
  }

  # ADDING EXTRA STUFF TO THE PLOTS

  df <- data.frame(x=coords, xend=total_length*(s-min(s))/(max(s)-min(s)), y=0, yend=min_height)
  # Control segment between diagram to exon
  # if (! is.null( exons_table) ){
  # plots[[length(plots)]] <- plots[[length(plots)]] + geom_segment(data=df, aes(x=x,y=y,xend=xend,yend=yend),alpha=0.2, lty=2)
  # }
  
  # adding the exons!
  
  if (!is.null(exons_table)) {
    exons_chr <- exons_table[exons_table$chr==intron_meta$chr[1],] # subset by chr
    exons_here <- exons_chr[ ( min(s) <= exons_chr$start & exons_chr$start <= max(s) ) | ( min(s) <= exons_chr$end & exons_chr$end <= max(s) ), ] # find exons
    
    #return(list(intron_meta, exons_here))
    # remove exons that don't actually start or end with a junction
    # and any repeated exons
    exons_here <-  unique( exons_here[ 
                      exons_here$end %in% intron_meta$start | 
                        exons_here$start %in% intron_meta$end, ]
    )
    # remove any exons larger than 500bp
    exons_here <- exons_here[ exons_here$end - exons_here$start <= 500 |
                                exons_here$end == min(intron_meta$start) |
                                exons_here$start == max(intron_meta$end),]
    # FOR EXON DIAGNOSTICS
    # exon_intervals <- Intervals(exons_here[,2:3])
    # plot(exon_intervals, use_points = FALSE)
    # return(exons_here)
    # 
    exons_here$gene_name=factor(exons_here$gene_name)
    
    print(exons_here)
   
    # gene_heights=min_height - ((1:length(levels(exons_here$gene_name)))-1.0) * abs(min_height) *  .15 # 0.05
    # heights=gene_heights[ as.numeric(exons_here$gene_name)] # .15
    # df <- data.frame(
    #   x=total_length*(exons_here$start-min(s))/(max(s)-min(s)), 
    #   xend=total_length*(exons_here$end-min(s))/(max(s)-min(s)), 
    #   y=heights, 
    #   yend=heights)

    # if (nrow(exons_here)>0) {
    #   plots[[length(plots)]] = plots[[length(plots)]] + 
    #     geom_segment(data=df, aes(x=x,y=y,xend=xend,yend=yend), alpha=1, size=2.5)  + 
    #     geom_hline(yintercept=min_height,alpha=1,size=0.7) + 
    #     geom_text(data=data.frame(x=my_xlim[1], y=gene_heights+1.2, label=levels(exons_here$gene_name)),aes(x,y,label=label))
    # }
    if (nrow( exons_here) > 0) {
      n_genes <- seq(1, length( levels(exons_here$gene_name) ) )
      gene_name_df <- data.frame( x= 0.2*total_length,#(n_genes * total_length) / (max(n_genes) + 1), 
                                  y=YLIMN, 
                                  label=rev(levels(exons_here$gene_name)))
      # count the occurences of each gene's exons and order by it
      gene_name_df$N <- table(exons_here$gene_name)[gene_name_df$label]
      gene_name_df <- gene_name_df[ order(gene_name_df$N, decreasing = TRUE),]
      
      #print(gene_name_df)
      
      # sort out colour pallette - most common gene goes first!
      cbbPalette <- cbbPalette[1:length(gene_name_df$label)]
      names(cbbPalette) <- gene_name_df$label
      # add junction colours 
      cbbPalette <- c(cbbPalette, junction_colour, cryptic_colour)
      names(cbbPalette)[(length(cbbPalette)-1):length(cbbPalette)] <- c("annotated","cryptic")
      #print(cbbPalette)
      plots[[1]] <- plots[[1]] + 
         guides(colour=FALSE) # don't show legend in top plot
        # geom_text( data = gene_name_df,
        #            fontface = "bold.italic",
        #            size = 6,
        #            nudge_x = 2,
        #            nudge_y = 1,
        #            aes(x,y,label=label,colour = label)
        # )
      plots[[2]] <- plots[[2]] +
        theme(legend.position="bottom", legend.justification = 'right') 

    }
    
    # pick principal gene - the one that contributes the most exons to the cluster will be at the top of the df
    principal_gene <- gene_name_df$label[1]
    
    #print(exons_here)
    gene_strand <- unique(exons_here[exons_here$gene_name == principal_gene,]$strand)
    #print(gene_strand)
    if( length(gene_strand)  > 1 ){
      gene_strand <- NA
    }
    if( !is.na(gene_strand) ){
      if( gene_strand == "+" ){
        strand_pos <- "last"
      }
      if( gene_strand == "-" ){
        strand_pos <- "first"
      }
    }else{
      strand_pos <- NULL
    }
    
    #print(strand_pos)
    
    invert_mapping=function(pos){
      if (pos %in% s) coords[as.character(pos)] else
        if (pos < min(s)) my_xlim[1] else
          if (pos > max(s)) my_xlim[2] else {
            w=which( pos < s[2:length(s)] & pos > s[1:(length(s)-1)] )
            stopifnot(length(w)==1)
            coords[w] + (coords[w+1]-coords[w])*(pos - s[w])/(s[w+1]-s[w])
          }
      }
      if (nrow(exons_here)>0) {
      exon_df <- data.frame( x=sapply(exons_here$start,invert_mapping), 
                        xend=sapply(exons_here$end,invert_mapping), 
                        y=0, 
                        yend=0,
                        label = exons_here$gene_name)
      
      # alter exon sizes to conform to a minimum exon length
      exon_df[ (exon_df$xend - exon_df$x) < min_exon_length, ]$xend <- exon_df[ (exon_df$xend - exon_df$x) < min_exon_length, ]$x + min_exon_length

      #print(exon_df)
      
      # find the negative space between exons
      # starts <- df$x[2:length(df$x)]
      # ends <- df$xend[1:length(df$xend) - 1]
      # spaces <- data.frame(x = ends, xend = starts)
      # print(spaces)
      #print(total_length)
      # print(nrow(df))
      # remove any exons that fall outside of the plotting space
      #df <- df[ df$x >= 0 & df$xend <= my_xlim[2], ]
      # print(df)
      # print(nrow(df))
      # for strand arrows
      if( !is.null(strand_pos)){
          # strand_df <- data.frame( x = c( -.05 * total_length,
          #                                 0.3 * total_length,
          #                                 0.3 * total_length,
          #                                 0.55 * total_length,
          #                                 0.55 * total_length,
          #                                 0.8 * total_length,
          #                                 0.8 * total_length,
          #                                 1.05 * total_length),
          #                          group = c( rep(1,2), rep(2,2), rep(3,2), rep(4,2) ),
          #                          y = 0)
        exon_intervals <- Intervals( matrix(data = c(exon_df$x, exon_df$xend), ncol = 2) )
        exon_intervals <- interval_union( exon_intervals )
        intron_intervals <- interval_complement( exon_intervals )
        intron_intervals <- intron_intervals[ 2:(nrow(intron_intervals)-1),]
        # strand arrows should be placed between exons
        #plot(intron_intervals, use_points=FALSE)
        strand_df <- as.data.frame(intron_intervals)
        # remove strand arrows where the introns are too small
        strand_df <- strand_df[ (strand_df$V2 - strand_df$V1) > 0.025 * total_length,]
        strand_df$midpoint <- strand_df$V1 + (strand_df$V2 - strand_df$V1) / 2
        
        group <- c()
        for(i in 1:nrow(strand_df)){
          group <- c(group, rep(i,2))
        }
        
        if( gene_strand == "+"){
          strand_df <- data.frame( x = c(rbind(strand_df$V1, strand_df$midpoint)),
                                   group = group, 
                                   y = 0) 
          
        }
        if( gene_strand == "-"){
          strand_df <- data.frame( x = c(rbind(strand_df$midpoint, strand_df$V2)),
                                   group = group,
                                   y = 0)
        }
        
        print(strand_df)
        
      # add strand arrows
        for (i in 1:length(plots) ){
          plots[[i]] <- plots[[i]] +
            geom_line( data = strand_df, 
                       aes( x = x, y = y, group = group ), colour = "black", size=1, 
                       arrow = arrow(ends = strand_pos, type = "open", angle = 30, length = unit(0.1, units = "inches" ))) 
        }
    }
    # add exons
    for (i in 1:length(plots) ){ 
      plots[[i]] <- plots[[i]] +
        scale_colour_manual("", values = cbbPalette ) +
        geom_segment(data=exon_df, aes(x=x,y=y,xend=xend,yend=yend, colour = label), alpha=1, size=6) +
        geom_segment( data = exon_df, aes(x = x - 0.01, xend = x, y = y, yend = yend), colour = "white", size = 6, alpha = 1) +
        geom_segment( data = exon_df, aes(x = xend, xend=xend + 0.01, y = y, yend = yend), colour = "white", size = 6, alpha = 1)
    }
    }
  }

  
  if (!is.na(main_title)) plots[[1]] = plots[[1]] + ggtitle(main_title)
  # arrange plots
  #do.call( gridExtra::grid.arrange, c(plots, list(ncol=1)))
  # only if no SNP is present
  if( is.na(snp_pos) ){
    gridExtra::grid.arrange( plots[[1]], plots[[2]], ncol =1)
  }
  #if (!is.null(exons_table))
  #exons_here

  }
  
