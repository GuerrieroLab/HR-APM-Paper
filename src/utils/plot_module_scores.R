my_plot_d1 <- function(data, x1){
  
  exp.pos <- (data %>% filter(mod %in% c("HR+/HER2-","HR+/HER2+")))[[x1]]
  exp.neg <- (data %>% filter(mod %in% c("HR-/HER2-","HR-/HER2+")))[[x1]]
  
  predictions <- c(exp.pos,exp.neg)
  labels <- rep(c(1,-1),c(length(exp.pos),length(exp.neg)))
  if(x1 %in% c("ERS","ESR1")){
    labels <- labels
  }else{
    labels <- -labels
  }
  pred <- prediction(predictions,labels)
  
  perf <- performance(pred, "tpr", "fpr")
  perf
  
  th <- threshold1(predictions,labels)
  
  p <- ggplot(data,aes(x=mod,y=!!sym(x1),fill=is.hr)) +
    geom_violin(position="dodge") +
    scale_fill_manual(values=brewer.pal(4,"YlOrBr")[c(3,1)])+
    geom_jitter(width=.2,size=.5) +
    theme_bw() +
    xlab("") +
    ylab("") +
    ggtitle(paste0(x1," activity")) +
    theme(plot.title=element_text(hjust=0.5),
          text = element_text(size = 20),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position="none")  
  
  p <- p + geom_hline(yintercept=th,color=1,linetype="dashed")
  
  return(p)  
}

my_plot_d2 <- function(data, x1="ERS", x2){
  
  exp.pos <- (data %>% filter(mod %in% c("HR+/HER2-","HR+/HER2+")))[[x1]]
  exp.neg <- (data %>% filter(mod %in% c("HR-/HER2-","HR-/HER2+")))[[x1]]
  
  predictions <- c(exp.pos,exp.neg)
  labels <- rep(c(1,-1),c(length(exp.pos),length(exp.neg)))
  if(x1 %in% c("ERS","ESR1")){
    labels <- labels
  }else{
    labels <- -labels
  }
  
  dens <- densCols(data[[x1]], data[[x2]], colramp=colorRampPalette(blues9[-(1:3)]))
  
  data$dens <- dens
  
  p <- ggplot(data,aes(x=!!sym(x1),y=!!sym(x2))) +
    geom_point(color=dens) +
    # scale_color_gradientn(colours=rainbow(100)) +
    theme_bw() +
    xlab(x1) +
    ylab(x2) +
    ggtitle(paste0(x1," vs ",x2)) +
    geom_smooth(color=2,span=.8)+
    theme(legend.position="none")
  return(p)  
}


my_plot_d3 <- function(data, x.th,x1="ERS", x2,draw.th=TRUE,split.lm=TRUE){
  
  if(missing(x.th)){
    x.th <- x1
  }
  exp.pos <- (data %>% filter(mod %in% c("HR+/HER2-","HR+/HER2+")))[[x.th]]
  exp.neg <- (data %>% filter(mod %in% c("HR-/HER2-","HR-/HER2+")))[[x.th]]
  
  predictions <- c(exp.pos,exp.neg)
  labels <- rep(c(1,-1),c(length(exp.pos),length(exp.neg)))
  if(x.th %in% c("ERS","ESR1")){
    labels <- labels
  }else{
    labels <- -labels
  }
  th <- threshold1(predictions,labels)
  
  data <- data %>%
    mutate(is.pos = !!sym(x.th) > th)
  
  ##
  dens <- densCols(data[[x1]], data[[x2]], colramp=colorRampPalette(blues9[-(1:3)]))
  
  data$dens <- dens
  
  p <- ggplot(data,aes(x=!!sym(x1),y=!!sym(x2))) +
    geom_point(color=dens) +
    # scale_color_gradientn(colours=rainbow(100)) +
    theme_bw() +
    xlab(x1) +
    ylab(x2) +
    ggtitle(paste0(x1," vs ",x2)) +
    # geom_smooth(aes(color = is.pos), method="lm",span=.8)+
    theme(legend.position="none")
  if(draw.th){
    p <- p + geom_vline(xintercept = th,linetype=2)
  }
  if(split.lm){
    p <- p +
      stat_smooth(method=function(formula,data,weights=weight) rlm(formula,
                                                                   data=data,
                                                                   weights=weight,
                                                                   method="MM"),
                  aes(color=is.pos))
  }else{
    p <- p +
      stat_smooth(method=function(formula,data,weights=weight) rlm(formula,
                                                                   data=data,
                                                                   weights=weight,
                                                                   method="MM"),
                  color = 1)
  }
  
  return(p)  
}


my_plot_d2_th <- function(data, x1="ERS", x2){
  
  exp.pos <- (data %>% filter(mod %in% c("HR+/HER2-","HR+/HER2+")))[[x1]]
  exp.neg <- (data %>% filter(mod %in% c("HR-/HER2-","HR-/HER2+")))[[x1]]
  
  predictions <- c(exp.pos,exp.neg)
  labels <- rep(c(1,-1),c(length(exp.pos),length(exp.neg)))
  if(x1 %in% c("ERS","ESR1")){
    labels <- labels
  }else{
    labels <- -labels
  }
  pred <- prediction(predictions,labels)
  
  perf <- performance(pred, "tpr", "fpr")
  
  th <- threshold1(predictions,labels)
  
  dens <- densCols(data[[x1]], data[[x2]], colramp=colorRampPalette(blues9[-(1:3)]))
  
  data$dens <- dens
  
  p <- ggplot(data,aes(x=!!sym(x1),y=!!sym(x2))) +
    geom_point(color=dens) +
    # scale_color_gradientn(colours=rainbow(100)) +
    theme_bw() +
    xlab(x1) +
    ylab(x2) +
    ggtitle(paste0(x1," vs ",x2)) +
    geom_smooth(color=2,span=.8)+
    geom_vline(xintercept=th,color=1,linetype="dashed") +
    theme(legend.position="none")
  return(p)  
}


my_plot_d3_th <- function(data, x.th,x1="ERS", x2,draw.th=TRUE,split.lm=TRUE){
  ## x1
  th <- function(x){
    exp.pos <- (data %>% filter(mod %in% c("HR+/HER2-","HR+/HER2+")))[[x]]
    exp.neg <- (data %>% filter(mod %in% c("HR-/HER2-","HR-/HER2+")))[[x]]
    
    predictions <- c(exp.pos,exp.neg)
    labels <- rep(c(1,-1),c(length(exp.pos),length(exp.neg)))
    
    if(x %in% c("ERS","ESR1")){
      labels <- labels
    }else{
      labels <- -labels
    }
    
    pred <- prediction(predictions,labels)
    
    perf <- performance(pred, "tpr", "fpr")
    
    th <- threshold1(predictions,labels)
    return(th)
  }
  
  th1 <- th(x1)
  th2 <- th(x2)  

  if(missing(x.th)){
    x.th <- x1
  }
  
  exp.pos <- (data %>% filter(mod %in% c("HR+/HER2-","HR+/HER2+")))[[x.th]]
  exp.neg <- (data %>% filter(mod %in% c("HR-/HER2-","HR-/HER2+")))[[x.th]]
  
  predictions <- c(exp.pos,exp.neg)
  labels <- rep(c(1,-1),c(length(exp.pos),length(exp.neg)))
  if(x.th %in% c("ERS","ESR1")){
    labels <- labels
  }else{
    labels <- -labels
  }
  th <- threshold1(predictions,labels)
  
  data <- data %>%
    mutate(is.pos = !!sym(x.th) > th)
  
  ##
  dens <- densCols(data[[x1]], data[[x2]], colramp=colorRampPalette(blues9[-(1:3)]))
  
  data$dens <- dens
  
  p <- ggplot(data,aes(x=!!sym(x1),y=!!sym(x2))) +
    geom_point(color=dens) +
    # scale_color_gradientn(colours=rainbow(100)) +
    theme_bw() +
    xlab(x1) +
    ylab(x2) +
    ggtitle(paste0(x1," vs ",x2)) +
    theme(legend.position="none")
  if(draw.th){
    p <- p + 
      geom_vline(xintercept = th1,linetype=2) +
      geom_hline(yintercept = th2,linetype=2)
  }
  if(split.lm){
    p <- p +
      stat_smooth(method=function(formula,data,weights=weight) rlm(formula,
                                                                   data=data,
                                                                   weights=weight,
                                                                   method="MM"),
                  aes(color=is.pos))
  }else{
    p <- p +
      stat_smooth(method=function(formula,data,weights=weight) rlm(formula,
                                                                   data=data,
                                                                   weights=weight,
                                                                   method="MM"),
                  color = 1)
  }
  
  return(p)  
}

my_plot_d2_th2 <- function(data, x1="ERS", x2){
  
  ## x1
  th <- function(x){
    exp.pos <- (data %>% filter(mod %in% c("HR+/HER2-","HR+/HER2+")))[[x]]
    exp.neg <- (data %>% filter(mod %in% c("HR-/HER2-","HR-/HER2+")))[[x]]
    
    predictions <- c(exp.pos,exp.neg)
    labels <- rep(c(1,-1),c(length(exp.pos),length(exp.neg)))
    
    if(x %in% c("ERS","ESR1")){
      labels <- labels
    }else{
      labels <- -labels
    }
    
    pred <- prediction(predictions,labels)
    
    perf <- performance(pred, "tpr", "fpr")
    
    th <- threshold1(predictions,labels)
    return(th)
  }
  
  th1 <- th(x1)
  th2 <- th(x2)
  
  dens <- densCols(data[[x1]], data[[x2]], colramp=colorRampPalette(blues9[-(1:3)]))
  
  data$dens <- dens
  
  p <- ggplot(data,aes(x=!!sym(x1),y=!!sym(x2))) +
    geom_point(color=dens) +
    # scale_color_gradientn(colours=rainbow(100)) +
    theme_bw() +
    xlab(x1) +
    ylab(x2) +
    ggtitle(paste0(x1," vs ",x2)) +
    geom_smooth(color=2,span=.8)+
    geom_vline(xintercept=th1,color=1,linetype="dashed") +
    geom_hline(yintercept=th2,color=1,linetype="dashed") +    
    theme(legend.position="none")
  return(p)  
}

###
th <- function(df,x){
  exp.pos <- (df %>% filter(mod %in% c("HR+/HER2-","HR+/HER2+")))[[x]]
  exp.neg <- (df %>% filter(mod %in% c("HR-/HER2-","HR-/HER2+")))[[x]]
  
  predictions <- c(exp.pos,exp.neg)
  labels <- rep(c(1,-1),c(length(exp.pos),length(exp.neg)))
  
  if(x %in% c("ERS","ESR1")){
    labels <- labels
  }else{
    labels <- -labels
  }
  
  pred <- prediction(predictions,labels)
  
  perf <- performance(pred, "tpr", "fpr")
  
  th <- threshold1(predictions,labels)
  return(th)
}

###
my_plot_d2_col <- function(data, x1="ERS", x2, x.col){
  
  th1 <- th(df=data,x1)
  th2 <- th(df=data,x2)
  
  # cols <- rev(colorRampPalette(brewer.pal(11,"RdBu"))(100))
  cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(100)
  
  p <- ggplot(data,aes(x=!!sym(x1),y=!!sym(x2))) +
    geom_point(aes(color=!!sym(x.col))) +
    scale_color_gradientn(colors=cols) + 
    theme_bw() +
    xlab(x1) +
    ylab(x2) +
    ggtitle(paste0(x1," vs ",x2)) +
    geom_vline(xintercept=th1,color=1,linetype="dashed") +
    geom_hline(yintercept=th2,color=1,linetype="dashed")
  theme(legend.position="none")
  return(p)  
}