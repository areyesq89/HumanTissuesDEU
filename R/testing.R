
estimateOneDispersionTRT <- function( modelMatrix, countVector,
                                     sizeFactors, dispInitialGuess = .5 ) {
    fit1 <- glmnb.fit( modelMatrix, countVector,
                      dispersion=dispInitialGuess, offset=log(sizeFactors) )
    exp( optimize( function( logalpha ){
                      -profileLogLikelihood(
                          exp(logalpha),
                          modelMatrix,
                          countVector,
                          fitted.values(fit1) ) },
                    log( c( 1e-5, 1e3 ) ) )$minimum )
}

profileLogLikelihood <- function( disp, mm, y, muhat ){
    if(length(disp) != length(y)){
        disp <- rep(disp, length(y))
    }
    ll <- sum( sapply( seq(along=y), function(i){
                          dnbinom( y[i], mu=muhat[i], size=1/disp[i], log=TRUE )
                      } ) )
    z <- log(muhat) + ( y - muhat ) / muhat
    v0 <- muhat + disp * muhat^2
    w <- 1 / ( ( 1 / muhat )^2 * v0 )
    qrres <- qr( mm*sqrt(w) )
    cr <- sum( log( abs( diag( qrres$qr )[ seq_len(qrres$rank) ] ) ) )
    ll - cr
}

#' @title Estimate dispersions
#' @description
#' This function estimates dispersions
#' for each exonic region of a \code{DEXSeqDataSet} object,
#' as in the \code{DEXSeq} method.
#' 
#' @param dxd A \code{DEXSeqDataSet} object
#' @param formula A formula to estimate dispersions
#' @param bjp A BiocParam instance
#'
#' @return A vector of relative exon usage coefficients
#' @export 
estimateTestDispersionsParallel <- function(dxd, formula, bjp){
    toSplit <- sort(rep(1:bpnworkers(bjp), length.out=nrow(dxd)))
    spMat <- split( as.data.frame(counts(dxd)), toSplit )
    mm <- DEXSeq:::rmDepCols( model.matrix(formula, as.data.frame(colData(dxd))))
    sizeFactors <- colData(dxd)$sizeFactor
    dispsAll <- bplapply( spMat,
        function(x, testableVectorZ, mmZ, sizeFactorsZ){
            disps <- sapply( rownames(x), function(ex){
                cat(sprintf("fitting %s\n", ex))
                if( ! testableVectorZ[ex]){
                    return(NA)
                }
                counts <- as.numeric(x[ex,])
                dp <- try( estimateOneDispersionTRT( mmZ, counts, sizeFactorsZ) )
                if( !inherits( dp, "try-error" ) ){
                    dp
                }else{
                    NA
                }
            } )
            names( disps ) <- rownames(x)
            disps
        }, testableVectorZ=testableVector, mmZ=mm,
        sizeFactorsZ=sizeFactors, BPPARAM=bjp )
    dispsAll <- unlist(dispsAll)
    dispsAll
}

#' @title Testing for differential exon usage
#' @description
#' This function tests for differential exon usage
#' for each exonic region of a \code{DEXSeqDataSet} object.
#' 
#' @param dxd A \code{DEXSeqDataSet} object
#' @param formula A formula to estimate dispersions
#' @param bjp A BiocParam instance
#'
#' @return A vector of p-values
#' @export 
testForDEUParallel <- function( dxd, null, full, rawDisps, testableVector, bjp ){
    modelFrame <- as.data.frame(colData(dxd))
    countMatrix <- counts(dxd)
    mmNull <- DEXSeq:::rmDepCols(model.matrix(null, modelFrame))
    mmFull <- DEXSeq:::rmDepCols(model.matrix(full, modelFrame))
    testableVector <- testableVector
    toSplit <- sort(rep(seq_len(bpnworkers(bjp)), length.out=nrow(countMatrix)))
    splitMat <- split( as.data.frame(countMatrix), toSplit )
    cat( sprintf("testing using %s cores\n", bpnworkers(bjp) ) )
    pvalAll <- bplapply( splitMat,
        function( mat, modelFrameZ, rawDispsZ,
                 mmNullZ, mmFullZ, testableVectorZ ){
            pvals <- sapply( rownames(mat), function(ex){
                cat(sprintf("testing %s\n", ex))
                x <- which(rownames(mat) %in% ex)
                xi <- which(rownames(rawDispsZ) %in% ex)
                if( !testableVectorZ[xi] ){
                    return(NA)
                }
                disps <- ifelse( modelFrameZ$exon=="this",
                    rawDispsZ[xi,"dispThis"],
                    rawDispsZ[xi,"dispOthers"] )
                counts <- as.numeric(mat[x,])
                sf <- modelFrameZ$sizeFactor
                options(warn=2)
                fitNull <- try( glmnb.fit(mmNullZ, counts, dispersion=disps, offset=sf) )
                fitFull <- try( glmnb.fit(mmFullZ, counts, dispersion=disps, offset=sf) )
                options(warn=0)
                if( inherits(fitNull, "try-error") | inherits( fitFull, "try-error") ){
                    return(NA)
                }
                pval <- 1 - pchisq( deviance(fitNull) - deviance(fitFull), df=ncol(mmFullZ) - ncol(mmNullZ))
                return(pval)
            } )
            names(pvals) <- rownames(mat)
            pvals
        }, modelFrameZ=modelFrame,
        rawDispsZ=rawDisps, mmNullZ=mmNull,
        mmFullZ=mmFull, testableVectorZ=testableVector,
        BPPARAM=bjp )
    unlist(pvalAll)
}
