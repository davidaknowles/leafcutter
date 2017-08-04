
#' Make gene level plot
#'
#' Take a gene name or ID, search the counts matrix for all clusters that fall within that gene's boundaries (found from the exon table). 
#' @import ggplot2
#' @export
make_gene_plot <- function(gene_name,
                           clusterID=NULL,
                           cluster_list=NULL,
                           introns=NULL,
                           introns_to_plot=NULL,
                           counts,
                           exons_table=NULL, 
                           len=500,
                           min_exon_length=0.5,
                           main_title=NA,
                           snp_pos=NA,
                           summary_func=colSums,
                           legend_title="Mean counts"){
  
  stopifnot( !is.null(exons_table) )
  stopifnot( !is.null(introns_to_plot))
  # only get the exons of the gene of interest - no overlapping genes!
  exons <- exons_table[ exons_table$gene_name == gene_name, ]
  
  gene_start <- min(exons$start)
  gene_end <- max(exons$end)
  
  stopifnot( length( unique(exons$chr) ) == 1 )
  
  # subset out the introns to plot based on whether they fall within the coordinates

  clusters <- introns_to_plot[ introns_to_plot$start >= gene_start & introns_to_plot$end <= gene_end & introns_to_plot$chr == unique(exons$chr) , ]
  cluster_ids <- introns_to_plot[ introns_to_plot$start >= gene_start & introns_to_plot$end <= gene_end& introns_to_plot$chr == unique(exons$chr) , ]$clu 
  
  #return(list(clusters,cluster_ids))
  
  stopifnot( nrow( clusters) > 0 )
  
  if( !is.null(clusterID)){
    stopifnot( clusterID %in% cluster_ids)
  }
  
  # if you want to represent the counts
  #y <- t(counts[grepl( paste(cluster_ids, collapse = "|"), introns_to_plot$clu), ])

  # remove duplicates
  exons <- exons[ !duplicated(paste( exons$start, exons$end)) ,]
  # remove exons over 500bp in length
  exons <- exons[ exons$end - exons$start <= 500,]
  # sort by end
  exons <- exons[ order(exons$end),]
  # for each cluster, remove any exon that fall with the cluster but are not connected by any junctions
  toDiscard <- c()
  for( clu in unique(cluster_ids) ){
    print(clu)
    clusterMin <- min(clusters[clusters$clu == clu,]$start)
    clusterMax <- max(clusters[clusters$clu == clu,]$end)
    clusterExons <- exons[ exons$end >= clusterMin & exons$start <= clusterMax,]
    notConnected <- clusterExons[ !(clusterExons$end %in% clusters[clusters$clu == clu,]$start |
                                 clusterExons$start %in% clusters[clusters$clu == clu,]$end), ]
    toDiscard <- c(toDiscard, paste(notConnected$start, notConnected$end))
    
  }
  to_keep=!(paste(exons$start, exons$end) %in% toDiscard)
  if (sum(to_keep) > 2)
    exons <- exons[ to_keep, ]
  
  # constitutive introns that connect each exon
  # NO LONGER USED - MAY COME BACK IN 
  # exon_junctions <- data.frame( start = c(NA, exons$end), end = c(exons$start, NA)  )
  # exon_junctions <- exon_junctions[ !is.na(exon_junctions$start) & !is.na(exon_junctions$end),]
  # # remove exon-exon junctions that are found in the clusters
  # exon_junctions <- exon_junctions[ !(exon_junctions$start %in% clusters$start) & !(exon_junctions$end %in% clusters$end), ]
  # 
  # even_junctions <- clusters[seq(2,nrow(clusters), 2),]
  # odd_junctions <- clusters[seq(1,nrow(clusters), 2),]

  # # plot the normal values - no scaling applied!

  
  # SCALE JUNCTIONS
  
  clusters$id <- "cluster"
  all_junctions <- clusters[, c("start","end", "clu")]
  
  # SCALE EXONS WITH SAME METHOD AS JUNCTIONS
  exons$clu <- row.names(exons)
  all_exons <- exons[, c("start", "end", "clu" )]
  #all_junctions <- rbind( all_junctions, all_exons)
  

  length_transform <- function(g) log(g+1)
  m <- reshape2::melt(all_junctions, id.vars = "clu") # melt to get list of all coordinates
  s <- unique(m$value) # get unique start and end values
  #if (!is.na(snp_pos)) s=c(s,snp_pos)
  s <- sort(s) 
  d <- s[2:length(s)] - s[1:length(s)-1] # get the difference between pairs of coordinates
  trans_d <- length_transform(d) # e.g. log(d+1), sqrt(d), d ^ .33 # apply a trasnformation function - doesn't work on negative numbers!
  coords <- c(0,cumsum(trans_d))
  names(coords) <- s
  total_length=sum(trans_d) # ==max(coords)
  # this is a problem - xlim is set to the total distance between the clusters - not the exons in the gene
  my_xlim=c(-.05*total_length,1.05*total_length)
  
  # SCALE EXON COORDINATES
  
  exons_here <- exons
  exons_here$gene_name=factor(exons_here$gene_name)

  # currently if a coordinate falls outside of the clusters then it is converted to the minimum or maximum
  # this removes exons that fall outside the clusters so fix this
  invert_mapping=function(pos, s, coords, xlim){
    if (pos %in% s) coords[as.character(pos)] else
      if (pos < min(s)) xlim[1] else
        if (pos > max(s)) xlim[2] else {
          # for which coordinate in s is pos closest to?
          w=which( pos < s[2:length(s)] & pos > s[1:(length(s)-1)] )
          stopifnot(length(w)==1)
          # to this number add the difference between that position and the next multiplied by 
          coords[w] + (coords[w+1]-coords[w])*(pos - s[w])/(s[w+1]-s[w])
        }
  }
  
  
  print(exons_here)
  exon_df <- data.frame( x=sapply(exons_here$start,FUN = function(X) invert_mapping(pos = X, s=s, coords = coords, xlim = my_xlim) ), 
                         xend=sapply(exons_here$end,FUN = function(X) invert_mapping(pos = X, s=s, coords = coords, xlim = my_xlim) ), 
                         y=numeric(nrow(exons_here)), 
                         yend=numeric(nrow(exons_here)),
                         label = exons_here$gene_name)
  
  
  
  # if exons fall upstream of the first cluster
  upstream_exons <- exons_here[ exons_here$start < min(s),]
  if( nrow(upstream_exons) > 0){
    upstream_s <- c( min(upstream_exons$start), max(upstream_exons$end))
    upstream_coords <- c( -length_transform( upstream_s[2] - upstream_s[1] ), 0) # fit all exons within this space
    names(upstream_coords) <- upstream_s
    upstream_df <- data.frame( x = sapply( upstream_exons$start, FUN = function(X) invert_mapping( pos = X, s = upstream_s, coords = upstream_coords, xlim = upstream_coords)),
                               xend = sapply( upstream_exons$end, FUN = function(X) invert_mapping( pos = X, s = upstream_s, coords = upstream_coords, xlim = upstream_coords)),
                               y = 0,
                               yend = 0,
                               label = upstream_exons$gene_name)
    # replace exon_df top and bottom with the new values
    exon_df[ 1:nrow(upstream_df),] <- upstream_df
  }
  
  # exons that are downstream of the clusters
  downstream_exons <- exons_here[ exons_here$end > max(s),]
  if( nrow(downstream_exons) > 0){
    downstream_s <- c( min(downstream_exons$start), max(downstream_exons$end))
    
    downstream_coords <- c( max(coords),  max(coords) + length_transform( downstream_s[2] - downstream_s[1] ) ) # fit all exons within this space
    names(downstream_coords) <- downstream_s
    downstream_df <- data.frame( x = sapply( downstream_exons$start, FUN = function(X) invert_mapping( pos = X, s = downstream_s, coords = downstream_coords, xlim = downstream_coords)),
                                 xend = sapply( downstream_exons$end, FUN = function(X) invert_mapping( pos = X, s = downstream_s, coords = downstream_coords, xlim = downstream_coords)),
                                 y = 0,
                                 yend = 0,
                                 label = downstream_exons$gene_name)                           
    exon_df[ ( ( nrow(exon_df) + 1) - nrow(downstream_df) ):nrow(exon_df) ,] <- downstream_df
  }

  # alter exon sizes to conform to a minimum exon length
  
  exon_df[ (exon_df$xend - exon_df$x) < min_exon_length, ]$xend <- exon_df[ (exon_df$xend - exon_df$x) < min_exon_length, ]$x + min_exon_length
  
  # create new x limit based on the upstream and downstream exons
  my_xlim <- c( min(exon_df$x), max(exon_df$xend))
  
  # create gene_length term
  gene_length <- max(exon_df$xend) - min(exon_df$x)
  
  # CREATE NEW THEME FOR PLOT
  
  new_theme_empty <- theme_bw(base_size = 15)
  new_theme_empty$line <- element_blank()
  new_theme_empty$rect <- element_blank()
  new_theme_empty$strip.text <- element_blank()
  new_theme_empty$axis.text <- element_blank()
  
 
  # HARD CODED PLOT SETTINGS 
  
  YLIMN=5 # 8
  YLIMP=-5  # -9
  # plot junctions as curves
  curv = 0.25 # normally 0.5
  curveMax = 1
  curveExponent = 2
  yOffset = 0
  # horizontal white line to clean up the edges of the curvess
  centreLineWidth = 3
  junction_colour <- "#d66464"
  # for gene name colouring (black by default, other colours for multiple genes)
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
 
  
  min_height=0
  max_height=-1

  # originally set as 0.65
  yFactor = 0.8
  
  
  
    
  
  # SCALE ODD JUNCTIONS
    
  allEdges=do.call(rbind,foreach (i=1:nrow(all_junctions)) %do% {
    #allEdges=do.call(rbind,foreach (i=1:2) %do% {
    if (i%%2==1) return(NULL)  # only care about the even numbered junctions?
    #if (intron_meta$counts[i]==0) return(NULL)
    start=coords[ as.character(all_junctions$start[i]) ]
    end=coords[ as.character(all_junctions$end[i]) ]
    l=end-start  # length of junction
    #edge = data.frame(bezier(x=c(start, (start + end)/2, end), y=c(0, (1+sqrt(l)) * ( (i %% 2)*2-1 ), 0),evaluation = len))
    h=(1+sqrt(l)) * ( (i %% 2)*2-1 )
    min_height=min(min_height,h)
    max_height=max(max_height,h)
    
    edge = data.frame(start, end)
    edge$startv <- all_junctions$start[i]
    edge$endv <- all_junctions$end[i]
    edge$start <- start  # why do this again? You already set them
    edge$end <- end
    #edge$log10counts=intron_meta$counts[i]+1
    
    # label each junction
    #edge$label=format(intron_meta$prop[i],digits=2, scientific=FALSE)
    #edge$label=format(intron_meta$counts[i], scientific=F)
    #edge$label=paste(format(intron_meta$counts[i], scientific=F),'\n(',format(intron_meta$prop[i],digits=2),')',sep='') 
    edge$clu<- all_junctions$clu[i]
    #edge$Curve <- j
    edge$Group <- i
    edge$xtext <-start+l/2
    edge$ytext <- -l^(yFactor)/2+0.5  # magic formula here
    #edge$SIZE <- intron_meta$prop[i]
    #edge$SIZE <- intron_meta$prop[i]+1 # make SIZE equal to the normalised count + 1
    #edge$SIZE <- as.factor(.bincode(intron_meta$prop[i]+0.1, breaks=seq(0,100,1)/100, include.lowest=TRUE))
    edge
  })
    
  # SCALE EVEN JUNCTIONS  
  
  allEdgesP=do.call(rbind,foreach (i=1:nrow(all_junctions)) %do% {
    if (i%%2==0) return(NULL)  # just the odd numbered junctions
    #if (intron_meta$counts[i]==0) return(NULL)
    start=coords[ as.character(all_junctions$start[i]) ]
    end=coords[ as.character(all_junctions$end[i]) ]
    l=end-start
    #edge = data.frame(bezier(x=c(start, (start + end)/2, end), y=c(0, (1+sqrt(l)) * ( (i %% 2)*2-1 ), 0),evaluation = len))
    h=(1+sqrt(l)) * ( (i %% 2)*2-1 )
    min_height=min(min_height,h)
    max_height=max(max_height,h)
    
    edge = data.frame(start, end)
    edge$startv <- all_junctions$start[i]
    edge$endv <- all_junctions$end[i]
    edge$start <- start
    edge$end <- end
    #edge$log10counts=intron_meta$counts[i]+1
    #edge$label=format(intron_meta$prop[i],digits=2, scientific=FALSE)
    #edge$label=format(intron_meta$prop[i],digits=2, scientific=T)
    #edge$label=paste(format(intron_meta$counts[i], scientific=F),'\n(',format(intron_meta$prop[i],digits=2),')',sep='')
    edge$clu<- all_junctions$clu[i]
    #edge$Curve <- j
    edge$Group <- i
    edge$xtext <-start+l/2
    edge$ytext <- l^(yFactor)/2-0.5
    #edge$SIZE <- intron_meta$prop[i]
    #edge$SIZE <- as.factor(.bincode(intron_meta$prop[i], breaks=seq(0,100,1)/100, include.lowest=TRUE))
    #edge$SIZE <- intron_meta$prop[i]+1
    edge
  })
  
  # if introns available then match deltaPSI 
  if( !is.null(introns) ){
    allEdges$deltaPSI <- introns$deltapsi[ match( paste(allEdges$clu, allEdges$startv, allEdges$endv ), paste( introns$clusterID, introns$start, introns$end) )]
    allEdgesP$deltaPSI <- introns$deltapsi[ match( paste(allEdgesP$clu, allEdgesP$startv, allEdgesP$endv ), paste( introns$clusterID, introns$start, introns$end) )]
    }
  
  
  # LABEL EACH CLUSTER
  
  junctions <- rbind( allEdges, allEdgesP)
  label_df <- as.data.frame( group_by(junctions, clu) %>% 
                               summarise( start = min(start), 
                                          end = max(end), 
                                          middle = start + ( (end - start) / 2) , 
                                          ytext = max(ytext) ) )
  
  
  
  if( !is.null(cluster_list) ){
    label_df$FDR <- cluster_list$FDR[ match( label_df$clu, cluster_list$clusterID)]
    label_df$FDR[ is.na(label_df$FDR)] <- "."
    label_df$label <- paste0( label_df$clu, "\n", label_df$FDR )
  }
    #print(label_df)
  
    # GENE HEIGHTS
    
    df <- data.frame(x=coords, xend=total_length*(s-min(s))/(max(s)-min(s)), y=0, yend=min_height)

    gene_heights=min_height - ((1:length(levels(exons_here$gene_name)))-1.0) * abs(min_height) *  .15 # 0.05
    heights=gene_heights[ as.numeric(exons_here$gene_name)] # .15
    df=data.frame(x=total_length*(exons_here$start-min(s))/(max(s)-min(s)), xend=total_length*(exons_here$end-min(s))/(max(s)-min(s)), y=heights, yend=heights)
    
    # GENE NAME
    n_genes <- seq(1, length( levels(exons_here$gene_name) ) )
    gene_name_df <- data.frame( x= total_length * 0.01, #(n_genes * total_length) / (max(n_genes) + 1), 
                                y=YLIMP - 0.1*YLIMP, 
                                label=levels(exons_here$gene_name))
    
    #print(exon_df)
    
    # STRANDING  
    
    gene_strand <- unique(exons_here$strand)
    if( length(gene_strand)  > 1 ){
      gene_strand <- NA
    }
    
    if( !is.na(gene_strand) ){
      strand_df <- exon_df[ exon_df$xend - exon_df$x >= 1, ]
      strand_df <- strand_df[ !duplicated( paste(strand_df$x, strand_df$xend) ), ]
      strand_df$id <- rownames(strand_df)
      # shorten the intervals slightly
      arrowFactor = 1
      strand_df$x = strand_df$x + arrowFactor 
      strand_df$xend <- strand_df$xend - arrowFactor
      # remove any arrows that are now too small
      strand_df <- strand_df[ strand_df$xend - strand_df$x >= 1, ]
      
      #print(strand_df)
      
      if( gene_strand == "+" ){
        strand_pos <- "last"
        # take first exon
        strand_df <- strand_df[ 1, ]
        # melt down into coordinates linked by id - the same exon
        strand_df <- melt( strand_df[, c("x","xend", "id")], "id")
        
      }
      if( gene_strand == "-" ){
        strand_pos <- "first"
        strand_df <- strand_df[ nrow(strand_df), ]
        # melt down into coordinates linked by id - the same exon
        strand_df <- melt( strand_df[, c("x","xend", "id")], "id")
        
      }
    }else{
      strand_pos <- NULL
    }
    
    # use intervals package to compute the correct placement of stranding arrows
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
    
    # STRAND ARROWS
    # only if exons are at least 1 unit in length
    
    #print(allEdgesP)
    #print(allEdges)
    
    # PLOTTING
    
    # colour the clusters with the lowest p-values
    # if( !is.null(clusterID) ){
    #   cluster_colours <- ifelse( allEdges$clu == clusterID, junction_colour, "gray" )
    #   cluster_coloursP <- ifelse( allEdgesP$clu == clusterID, junction_colour, "gray" )
    # }
    
    # colour only the clusters with significant p values
    # cluster_colours <- ifelse( label_df$FDR[ match(allEdges$clu, label_df$clu)] != ".", junction_colour, "gray" )
    # cluster_coloursP <- ifelse( label_df$FDR[ match(allEdges$clu, label_df$clu)] != ".", junction_colour, "gray" )
    
    #label_df$FDR[ match(allEdges$clu, label_df$clu)] != "."
    
    #print(exon_df)
    
    
    plots <- ggplot() + 
      # cluster junctions

      geom_curve(data=allEdgesP[ label_df$FDR[ match(allEdgesP$clu, label_df$clu)] != "." ,], aes(x = start, xend = end, y = 0, yend = 0, group = Group, color=clu, size = curveMax ),
                 curvature=curv,lineend="round", colour = junction_colour ) +
      geom_curve(data=allEdges[ label_df$FDR[ match(allEdges$clu, label_df$clu)] != "." ,], aes(x = start, xend = end, y = 0, yend = 0, group = Group, color=clu, size = curveMax),
                 curvature=-curv,lineend="round",  colour = junction_colour )
      
      if( nrow(  allEdgesP[ label_df$FDR[ match(allEdgesP$clu, label_df$clu)] == "." ,]) > 0  ){
      # if there are non-significant clusters, then colour  grey
      plots <- plots +
      geom_curve(data=allEdgesP[ label_df$FDR[ match(allEdgesP$clu, label_df$clu)] == "." ,], aes(x = start, xend = end, y = 0, yend = 0, group = Group, color=clu, size = curveMax ),
                 curvature=curv,lineend="round", colour = "gray" )
      }
      
      if( nrow(  allEdges[ label_df$FDR[ match(allEdges$clu, label_df$clu)] == "." ,]) > 0  ){
      plots <- plots +
      geom_curve(data=allEdges[ label_df$FDR[ match(allEdges$clu, label_df$clu)] == "." ,], aes(x = start, xend = end, y = 0, yend = 0, group = Group, color=clu, size = curveMax),
                 curvature=-curv,lineend="round",  colour = "gray" )
      }
      
      # to increase visual contrast convert all positive deltaPSI to 0.5 and all negative to -0.5
    
      if( !is.null(introns)){
        plots <- plots +
          geom_curve(data=allEdgesP[ label_df$FDR[ match(allEdgesP$clu, label_df$clu)] != "." ,], aes(x = start, xend = end, y = 0, yend = 0, group = Group, color=as.factor(deltaPSI > 0) , size = curveMax ),
                     curvature=curv,lineend="round") +
          geom_curve(data=allEdges[ label_df$FDR[ match(allEdges$clu, label_df$clu)] != "." ,], aes(x = start, xend = end, y = 0, yend = 0, group = Group, color=as.factor(deltaPSI > 0), size = curveMax),
                     curvature=-curv,lineend="round" ) +
          scale_colour_discrete("dPSI",labels = c("down","up") )
      }
      
     
      # geom_curve(data=allEdgesP, aes(x = start, xend = end, y = 0, yend = 0, group = Group, colour=cluster_coloursP, size = curveMax ),
      #            curvature=curv,lineend="round" ) +
      # geom_curve(data=allEdges,aes(x = start, xend = end, y = 0, yend = 0, group = Group, colour=cluster_colours, size = curveMax),
      #            curvature=-curv,lineend="round" ) +
      plots <- plots +      
      
      new_theme_empty + 
      # make the y axis label the group
      #ylab(paste0(tis," (n=",group_sample_size,")")) + 
      xlab("") + 
      ylab("") +
      xlim(my_xlim) +
      # try titling instead - why doesn't this work?
      #ggtitle(paste0(tis," (n=",group_sample_size,")" ) ) + 
      
      # horizontal line - smooth out the ends of the curves
      geom_hline(yintercept=0, size = centreLineWidth, colour = "white") +
      geom_hline(yintercept=0,alpha=.9, size=1) + 
      
      # label the junctions
      #geom_label(data=allEdgesP,aes(x=xtext,y=ytext,label=label), size = 5, label.size = 0 ) +  
      #geom_label(data=allEdges,aes(x=xtext,y=ytext,label=label), size= 5, label.size = 0 ) +
      
      ylim(YLIMN,YLIMP) + 
        
      scale_size_continuous(limits=c(0,10),guide='none') + 
      
      # is this used for anything? color is currently set to clu which doesn't change for each junction
      #scale_color_discrete( c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928'),guide='none') #+ scale_color_gradient2(low='white',mid='lightblue',high='purple4',guide="none")  #  + scale_color_discrete(guide="none")
      #scale_colour_manual( values = cbbPalette ) +
      #guides(colour=FALSE) + # don't show legend
# 
#       geom_text( data = gene_name_df, 
#                    fontface = "bold.italic",  
#                    size = 6, 
#                    nudge_x = 0, 
#                    nudge_y = 0, 
#                    aes(x,y,label=label ),colour = "black"
#         ) + 
        ggtitle( paste(gene_name_df$label, collapse = "+") ) +
        theme(plot.title = element_text(face="bold.italic", colour="black", size = 20)) +
      # EXONS
      geom_segment( data=exon_df, aes(x=x,y=y,xend=xend,yend=yend ), alpha=1, size=6, colour = 'black' )      
      
      #geom_segment( data = df, aes(x = x - 0.05, xend = x, y = y, yend = yend), colour = "white", size = 6, alpha = 1) +
      #geom_segment( data = df, aes(x = xend, xend=xend + 0.05, y = y, yend = yend), colour = "white", size = 6, alpha = 1)
    

    
      # ADD STRANDING ARROWS
    #print(paste("strand is: ", strand_pos) )
    print(strand_df)
    if( !is.null(strand_pos)){
        plots <- plots +
          geom_line( data = strand_df,
                     aes( x = x, y = y, group = group ),
                     colour = "black", size=1,
                     arrow = arrow(ends = strand_pos, type = "open", angle = 30, length = unit(0.1, units = "inches" ) )
                     )
      }
    # DEMARCATE JUNCTIONS AND EXONS
    
    plots <- plots +
      # # JUNCTIONS
      # geom_segment( data = allEdgesP, aes(x = end, xend=end + 0.05, y = 0, yend = 0), colour = "white", size = 6 ) +
      # geom_segment( data = allEdges, aes(x = end, xend=end + 0.05, y = 0, yend = 0), colour = "white", size = 6 ) +
      # geom_segment( data = allEdgesP, aes(x = start, xend=start - 0.05, y = 0, yend = 0), colour = "white", size = 6 ) +
      # geom_segment( data = allEdges, aes(x = start, xend=start - 0.05, y = 0, yend = 0), colour = "white", size = 6 ) +
      # EXONS
      geom_segment( data = exon_df, aes(x = xend, xend=xend + 0.05, y = 0, yend = 0), colour = "white", size = 6 ) +
      geom_segment( data = exon_df, aes(x = x, xend=x - 0.05, y = 0, yend = 0), colour = "white", size = 6 ) 
    
    
  
    
    # LABEL CLUSTERS
    label_max <- YLIMN - 0.1*YLIMN  
    plots <- plots +
      geom_text( data = label_df, aes( x = middle, y = label_max, label = label)) +
      theme(legend.position="bottom", legend.justification = 'right') 
    
    
    
   #print(label_df)
    
  #print(clusterID)  
    
    if( !is.null(clusterID) ){
      plots <- plots + 
        geom_label( data = label_df[ label_df$clu == clusterID,],  aes( x = middle, y = label_max, label = label), fontface = "bold", label.size = 0.5, label.r = unit(0.3,"lines"), label.padding = unit(0.5,"lines") )
      
    }
    print("gene plot:")    
    print(exon_df)
   plots
}
  #if (!is.na(main_title)) plots[[1]] = plots[[1]] + ggtitle(main_title)
  # arrange plots
  #do.call( gridExtra::grid.arrange, c(plots, list(ncol=1)))

  #if (!is.null(exons_table))
  #exons_here
  




