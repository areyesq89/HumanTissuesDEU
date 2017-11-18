
getIndexList <- function(dxd){
    modelFrame <- as.data.frame(colData(dxd))
    sampleIdx <- data.frame(
        sample = levels( modelFrame$sample ),
        col = seq_along( levels( modelFrame$sample ) ) )
    exonAvgIdx <- max(sampleIdx$col) + 1
    sexIdx <- exonAvgIdx + 1
    crossIdx <- expand.grid(
        tissue=levels(modelFrame$tissue),
        individual=levels(modelFrame$individual) )
    crossIdx$col = sexIdx + 1:nrow(crossIdx)
    list(
        sampleIdx=sampleIdx,
        exonAvgIdx=exonAvgIdx,
        sexIdx=sexIdx,
        crossIdx=crossIdx )
}

modelMatrixREUC <- function(dxd, indexList){
    modelFrame <- as.data.frame( colData(dxd) )
    mm <- matrix( 0, nrow = nrow(modelFrame), ncol = max(indexList[["crossIdx"]]$col) )
    colnames(mm) <- rep( NA, ncol(mm) )
    colnames(mm)[indexList[["sampleIdx"]]$col] <-
        paste( "sample", indexList[["sampleIdx"]]$sample, sep="_")
    colnames(mm)[indexList[["exonAvgIdx"]]] <- "exonAvg"
    colnames(mm)[indexList[["sexIdx"]]] <- "sex"
    colnames(mm)[indexList[["crossIdx"]]$col] <-
        sprintf( "tissue_%s:individual_%s", indexList[["crossIdx"]]$tissue, indexList[["crossIdx"]]$individual )
    for( i in 1:nrow(modelFrame) ) {
        mm[ i, indexList[["sampleIdx"]]$col[ indexList[["sampleIdx"]]$sample == modelFrame$sample[i] ] ] <- 1
        mm[ i, indexList[["exonAvgIdx"]] ] <- { if( modelFrame$exon[i] == "this" ) 1 else 0 }
        mm[ i, indexList[["sexIdx"]] ] <- { if( modelFrame$sex[i] == "male" && modelFrame$exon[i] == "this" ) .5 else -.5 }
        mm[ i, indexList[["crossIdx"]]$col[ indexList[["crossIdx"]]$tissue == modelFrame$tissue[i] & indexList[["crossIdx"]]$individual == modelFrame$individual[i] ] ] <- { if( modelFrame$exon[i] == "this" ) .5 else -.5 }
    }
    mm
}

ll <- function( muhat, y, disp ){
    sum( y * log( muhat ) - ( y + 1/disp ) * log( muhat + 1/disp ) )
}

dll <- function( muhat, y, disp ){
    y / muhat - ( y + 1/disp ) / ( muhat + 1/disp )
}

shrinkageFit <- function( mm, counts, dispersions,
                         sizeFactors, beta0, shrink, priorsd ) {
    ofit <- optim( beta0,
                  function(beta){
                      muhat <- exp( mm %*% beta ) * sizeFactors
                      -ll( muhat, counts, dispersions ) +
                          sum( beta[shrink]^2 ) / priorsd^2 / 2 } ,
                  function(beta){
                      muhat <- exp( mm %*% beta ) * sizeFactors
                      -t( dll( muhat, counts, dispersions ) * muhat ) %*% mm +
                          beta*shrink / priorsd^2 },
                  method="L-BFGS-B",
                  control = list( trace=0, maxit=5000, factr=1e5 ) )
    if( ofit$convergence != 0 )
        warning( "L-BFGS optimization did not converge." )
    ofit$par
}

llCR <- function( muhat, y, disp, mm ){
    ll <- sum( sapply( seq(along=y), function(i){
        dnbinom( y[i], mu=muhat[i], size=1/disp[i], log=TRUE ) } ) )
    z <- log(muhat) + ( y - muhat ) / muhat
    v0 <- muhat + disp * muhat^2
    w <- 1 / ( ( 1 / muhat )^2 * v0 )
    qrres <- qr( mm*sqrt(w) )
    cr <- sum( log( abs( diag( qrres$qr )[ seq_len(qrres$rank) ] ) ) )
    ll - cr
}

estimateDispForExon <- function( muhat, y, mm, isThis,
                                startDispThis=.1, startDispOthers=.1 ) {
    a <- optim( log( c( startDispThis, startDispOthers ) ),
               function( x ) {
                   disps <- ifelse( isThis, exp(x[1]), exp(x[2]) )
                   -llCR( muhat, y, disps, mm )
               } )
    names(a$par) <- c( "dispThis", "dispOthers" )
    exp(a$par)
}
