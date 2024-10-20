subtype_violin_within <- function(dfx,mod){
  dfx1 <- dfx %>% select(-2) %>% filter(type == mod) %>% unique()
  p0 <- ggplot(dfx %>% filter(type==mod),aes(x=subtype,y=sp,fill=is.hr)) +
    geom_hline(yintercept=0) + 
    geom_violin(position="dodge") +
    scale_fill_manual(values=brewer.pal(4,"YlOrBr")[c(3,1)])+
    theme_bw() +
    ggtitle(paste0("Within ",mod)) +
    geom_text(data=dfx1,aes(x=subtype,label=txt),y=1.1,size=4) + 
    scale_y_continuous(limits=c(-.7,1)) +
    ylim(c(-.7,1.2)) +
    xlab("") +
    ylab("SCC within module") +
    theme(legend.position="none",
          plot.title=element_text(hjust=0.5,size=12),
          text = element_text(size = 12),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))   
  return(p0)
}

subtype_violin_between <- function(dfx,mod){
  if(!mod %in% c("TNFA.NFKB", "IFN.I", "APM.TC")){
    stop("mod should be one of the immune module names")
  }
  mods <- paste0("IMM",1:3)
  mod1 <- mods[match(mod,c("TNFA.NFKB", "IFN.I", "APM.TC"))]
  dfx1 <- dfx %>% select(-2) %>% filter(type == paste0("ERS.",mod1)) %>% unique()
  
  p0 <- ggplot(dfx %>% filter(type==paste0("ERS.",mod1)),
               aes(x=subtype,y=sp,fill=is.hr)) +
    geom_hline(yintercept=0) + 
    geom_violin(position="dodge") +
    scale_fill_manual(values=brewer.pal(4,"YlOrBr")[c(3,1)])+
    theme_bw() +
    ggtitle(paste0("ERS vs ",mod)) +
    geom_text(data=dfx1,aes(x=subtype,label=txt),y=1.1,size=4) + 
    scale_y_continuous(limits=c(-.7,1)) +
    ylim(c(-.7,1.2)) +
    xlab("") +
    ylab("SCC between modules") +
    theme(legend.position="none",
          plot.title=element_text(hjust=0.5,size=12),
          text = element_text(size = 12),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))   
  return(p0)
}