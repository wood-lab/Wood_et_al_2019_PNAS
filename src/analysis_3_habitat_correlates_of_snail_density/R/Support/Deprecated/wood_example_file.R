panther.gg<-function(panther){
  require(ggplot2)
  p<-ggplot(panther,aes(x=TL_mm,y=philo_many_replaced,group=oldnew,shape=oldnew,col=oldnew,fill=oldnew))+

    scale_color_manual(values=c("darkblue","red"))+
    scale_fill_manual(values=c("darkblue","red"))+

    geom_point()+

    geom_smooth(aes(group=oldnew,color=oldnew),method="glm",family="poisson",se=TRUE)+

    xlab("total length of fish (mm)")+
    ylab("number of worms per fish")+

    theme(axis.title.y=element_text(size=20),axis.title.x=element_text(size=20),axis.text.y=element_text(size=10),axis.text.x=element_text(size=12),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),legend.position="none")

  return(p)
}
panther_plot<-panther.gg(panther)
panther_plot
