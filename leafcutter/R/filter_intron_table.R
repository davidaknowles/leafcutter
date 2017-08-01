
#' Filter intron table
#'
#' @import dplyr
#' @export 
filter_intron_table <- function(introns, clu, toSave=FALSE){
  d <- dplyr::filter(introns, clusterID == clu) %>% 
    dplyr::select( -clusterID, -gene, -ensemblID, -transcripts) %>% 
    arrange( desc(abs(deltapsi)))
  if( !toSave ){
   d <- rename(d, "ΔPSI" = deltapsi )
  }else{
    d <- rename(d, "dPSI" = deltapsi ) # fudge as grid arrange doesn't like greek letters
  }
  row.names(d) <- letters[1:nrow(d)] # letters is just a:z
  return(d)
}