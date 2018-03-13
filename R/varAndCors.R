#' @title Estimate coefficient of partial determination
#' @description
#' Estimate the coefficients of partial determination
#' 
#' @param crossCoefs A \code{DEXSeqDataSet} object
#' @param statistic A character vector with the name of the statistic (REUC or RSIC)
#' @param subset A character vector with the name of the subset.
#' @param BPPARAM BiocParallel instance
#'
#' @return A data frame with coefficients of partial determination
#' @export 

estimateR2 <- function( crossCoefs, statistic, subset, BPPARAM=SerialParam() ){         
    df <- data.frame(                                                                   
        individual=rep( rownames(crossCoefs[1,,]), ncol( crossCoefs[1,,] ) ),           
        tissue=rep( colnames(crossCoefs[1,,]), each=nrow( crossCoefs[1,,] ) ),          
        stringsAsFactors=FALSE )                                                        
    df$reuc <- NA                                                                       
    r2 <- bplapply( seq_along(dimnames(crossCoefs)[["exon"]]),                          
                   function(x, dfL){
                       df$reuc <- as.vector(crossCoefs[x,,])
                       if( any( is.na( df$reuc ) ) ){
                           return(rep(NA, 3))
                       }else{
                           anovaRes <- anova( lm( reuc ~ individual + tissue, df ) )
                           beta <- anovaRes[c("individual", "tissue", "Residuals"),"Sum Sq"]
                           return( beta )
                       }
                   }, dfL=df, BPPARAM=BPPARAM )
    res <- do.call(rbind, r2)                                                           
    colnames(res) <- c("individual", "tissue", "residuals")
    rownames(res) <- dimnames(crossCoefs)[["exon"]]
    res <- as.data.frame( res )
    res$statistic <- statistic
    res$subset <- subset
    res$exon <- rownames(res)                                                           
    res
}
