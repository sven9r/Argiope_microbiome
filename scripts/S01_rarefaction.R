#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr); library(ggplot2)
  library(vegan); library(ggpubr); library(readr); library(tibble)
})
source(here::here("scripts", "00_setup.R"))
raw <- load_argiope_raw()
bac <- raw$bacteria; coi <- raw$prey
get_otu <- function(x) grep("^OTU", names(x), value = TRUE)
mat_for <- function(data, otus, status = NULL) {
  if (!is.null(status)) data <- filter(data, Status == status)
  m <- data %>% select(ID, all_of(otus)) %>% mutate(across(all_of(otus), ~replace_na(as.numeric(.x), 0)))
  m <- m %>% filter(rowSums(across(all_of(otus))) > 1)
  mm <- as.matrix(select(m, all_of(otus))); storage.mode(mm) <- "numeric"
  rownames(mm) <- m$ID
  mm[, colSums(mm) > 0, drop = FALSE]
}
curves <- function(m, color, threshold, title, combined = FALSE, label_x, label_y) {
  grDevices::pdf(NULL)
  rc <- vegan::rarecurve(m, step = 100, label = FALSE)
  grDevices::dev.off()
  df <- bind_rows(lapply(seq_along(rc), function(i) tibble(ID=rownames(m)[i], reads=as.numeric(sub("^N", "", names(rc[[i]]))), richness=as.numeric(rc[[i]]))))
  if (combined) {
    status <- bind_rows(bac %>% select(ID,Status), coi %>% select(ID,Status)) %>% distinct(ID,.keep_all=TRUE)
    df <- left_join(df,status,by="ID")
    p <- ggplot(df,aes(reads,richness,group=ID,color=Status))+geom_line(linewidth=.6)+scale_color_manual(values=color)
  } else p <- ggplot(df,aes(reads,richness,group=ID))+geom_line(linewidth=.6,color=color)
  p + geom_vline(xintercept=threshold,color="darkgrey",linewidth=1)+theme_classic()+
    labs(title=title,x="Reads [n]",y="OTUs [n]",color="Age:")+theme(axis.text=element_text(size=11),axis.title=element_text(size=12),plot.title=element_text(size=13),legend.position=if(combined)c(.8,.8) else "none")+
    annotate("text",x=label_x,y=label_y,label=paste0("Cutoff: ",threshold," reads"),size=4)
}
bc <- mat_for(bac,get_otu(bac),"Contemporary"); pc <- mat_for(coi,get_otu(coi),"Contemporary")
ba <- mat_for(bac,get_otu(bac)); pa <- mat_for(coi,get_otu(coi))
p1 <- curves(bc,"#084594",1500,"Rarefaction curve of contemporary samples (16S)",label_x=7500,label_y=60)
p2 <- curves(pc,"#005A32",1693,"Rarefaction curve of contemporary samples (COI)",label_x=12500,label_y=8)
p3 <- curves(ba,c("Contemporary"="#2171B5","Museum"="#9ECAE1"),507,"Rarefaction curve of museum and contemporary samples (16S)",TRUE,7500,30)
p4 <- curves(pa,c("Contemporary"="#238B45","Museum"="#A1D99B"),1133,"Rarefaction curve of museum and contemporary samples (COI)",TRUE,15000,20)
plot <- ggarrange(p1,p2,p3,p4,nrow=2,ncol=2,labels=c("A","B","C","D"))
out <- here::here("figures", "supplement", "S01_rarefaction_curves.png")
ggsave(out,plot,width=14,height=8,dpi=300,bg="white")
write_csv(tibble(panel=c("A","B","C","D"),n_samples=c(nrow(bc),nrow(pc),nrow(ba),nrow(pa)),cutoff=c(1500,1693,507,1133)),here::here("output", "tables", "S01_rarefaction_panel_counts.csv"))
cat(out,"\n")
