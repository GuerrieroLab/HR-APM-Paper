tri.heatmap <- function(x,rm.diag=TRUE,mar=c(5,5,3,3),...){
  if(!isSymmetric(tcga.sp1)){
    stop("only symmetric matrix can be handled")
  }
  op <- par()
  par(mar=mar)
  
  n <- nrow(x)
  x[lower.tri(x)] <- NA
  if(rm.diag){
    m <- n-1
    idx <- seq(n)[-1]
    diag(x) <- NA
  }else{
    m <- n
    idx <- seq(n)
  }
  # cols <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
  cols <- colorRampPalette(c("navy", "white", "firebrick3"))(1024)
  image(seq(m),seq(m),x[seq(m),rev(seq(n)[idx])],col=cols,axes=F,xlab="",ylab="",zlim=c(-1,1),...)
  axis(1,at=seq(m),rownames(x)[seq(m)],las=2)
  axis(2,at=seq(m),rev(rownames(x)[idx]),las=1)
  for(i in seq(n)){
    for(j in seq(n)){
      text(i,n-j+1,round(x[i,j],2))
    }
  }
  on.exit(suppressWarnings(par(op)))
}