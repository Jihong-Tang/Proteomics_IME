library(tidyverse)

library(tidyverse)
pair_all <- readxl::read_xlsx("../Fig_codes/data/Fig5/0324_PA_Longitunidal_compare_AI_Model.xlsx")


plt_IR_diff <- pair_all[, c("PairID", "mean_ini2", "mean_rec2")] #%>% filter(group == "ProMES") #%>% filter(abs(diff) > 0.2)
#plt_IR_nodiff <- plt_IRcomapre %>% filter(abs(diff) <= 0.2)
long_data <- pivot_longer(plt_IR_diff, cols = -PairID, names_to = "IR", values_to = "value")

long_data$IR[long_data$IR == "mean_ini2"] <- "I"
long_data$IR[long_data$IR == "mean_rec2"] <- "R"

long_data$logvalue <- log(long_data$value * 100 +1)

Fig5B <- ggplot() + theme_classic()
Fig5B <- Fig5B + geom_boxplot(data = long_data, aes(x = IR, y = value),width = 0.4,size=0.8, alpha = 0.7,fill="transparent", outlier.colour = "white")
Fig5B
#Fig5B <- Fig5B + geom_line(data = plt_IR_nodiff, aes(x = IR, y = value,group=Patient_ID),size=0.25, color = "#dddddd", lty = "dashed")

#Fig5B <- Fig5B + geom_point(data = plt_IR_nodiff, aes(x = IR, y = value, color=nmf_subtype_I,group=Patient_ID),size=2, alpha = 0.7,shape=16) + 
#facet_wrap(factor(group, levels = c("ProMIX", "ProPPR", "ProMES", "ProNEU")) ~., ncol = 4)
Fig5B <- Fig5B + geom_line(data = long_data, aes(x = IR, y = value,group=PairID),size=0.25, color = "black",  lty = "dashed")
Fig5B <- Fig5B + geom_point(data = long_data, aes(x = IR, y = value, fill=IR,group=PairID),size=2, alpha = 0.7,stroke = .25,shape=21, color = "black")

Fig5B <- Fig5B + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                       plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                       legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                       axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                       axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))

Fig5B <- Fig5B + scale_y_continuous(expand=c(0,0),limits=c(-0.05,1),breaks = seq(-4,4,0.2)) +
  scale_x_discrete(limits = c("I","R"))
Fig5B<- Fig5B +ylab("NES") +xlab(NULL) 
Fig5B <- Fig5B +   scale_fill_manual(values = c("I" = "#009BFF", "R" = "#FF5D49")) 
#Fig5B <- Fig5B + scale_color_manual(name=NULL,values =  c("ProteoNMF1" = "#33a02c", "ProteoNMF2" = '#e31a1c', "ProteoNMF3" = '#ff7f00', "ProteoNMF4" = '#1f78b4', "Mixed" = "grey"))
Fig5B
figure_2<-rbind(ggplotGrob(Fig5B),size="last")
ggsave(file="./figures/Fig5/1127_boxplot_compare_ProSig_scores.pdf", plot=figure_2,bg = 'white', width =10, height = 7.5, units = 'cm', dpi = 600)
ggsave(file="./figures/Fig2/1115_boxplot_compare_ProSig_scores.pdf", plot=figure_2,bg = 'white', width =9, height = 15, units = 'cm', dpi = 600)


Fig5B <- ggplot() + theme_classic()
Fig5B <- Fig5B + geom_boxplot(data = long_data, aes(x = IR, y = logvalue),width = 0.4,size=0.8, alpha = 0.7,fill="transparent", outlier.colour = "white")
Fig5B
#Fig5B <- Fig5B + geom_line(data = plt_IR_nodiff, aes(x = IR, y = logvalue,group=Patient_ID),size=0.25, color = "#dddddd", lty = "dashed")

#Fig5B <- Fig5B + geom_point(data = plt_IR_nodiff, aes(x = IR, y = logvalue, color=nmf_subtype_I,group=Patient_ID),size=2, alpha = 0.7,shape=16) + 
#facet_wrap(factor(group, levels = c("ProMIX", "ProPPR", "ProMES", "ProNEU")) ~., ncol = 4)
Fig5B <- Fig5B + geom_line(data = long_data, aes(x = IR, y = logvalue,group=PairID),size=0.25, color = "black",  lty = "dashed")
Fig5B <- Fig5B + geom_point(data = long_data, aes(x = IR, y = logvalue, fill=IR,group=PairID),size=2, alpha = 0.7,stroke = .25,shape=21, color = "black")

Fig5B <- Fig5B + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                       plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                       legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
                       axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                       axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))

Fig5B <- Fig5B + scale_y_continuous(expand=c(0,0),limits=c(-0.05,1),breaks = seq(-4,4,0.2)) +
  scale_x_discrete(limits = c("I","R"))
Fig5B<- Fig5B +ylab("NES") +xlab(NULL) 
Fig5B <- Fig5B +   scale_fill_manual(values = c("I" = "#009BFF", "R" = "#FF5D49")) 
#Fig5B <- Fig5B + scale_color_manual(name=NULL,logvalues =  c("ProteoNMF1" = "#33a02c", "ProteoNMF2" = '#e31a1c', "ProteoNMF3" = '#ff7f00', "ProteoNMF4" = '#1f78b4', "Mixed" = "grey"))
Fig5B
figure_2<-rbind(ggplotGrob(Fig5B),size="last")
ggsave(file="./figures/Fig5/1127_boxplot_compare_ProSig_scores.pdf", plot=figure_2,bg = 'white', width =10, height = 7.5, units = 'cm', dpi = 600)

ggplot(plt_IR_diff) + 
  geom_point(aes(x = mean_ini2, y = mean_rec2))+
  theme_classic()



library(ggpubr)
long_data %>% 
  ggplot()+
  geom_boxplot(aes(x = IR, y =value),width = 0.5)+ 
  geom_line(aes(x = IR, y = value,group = PairID), size=1, color='gray', alpha=0.6)+ 
  geom_point(aes(x = IR, y = value,group=PairID),size=2,shape=16)+
  stat_compare_means(aes(x=IR, y = value), paired = T)+
  scale_fill_manual(values = c("#1f78b4", "#e31a1c"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "", y = "ProMES score", fill = "")

# 4D model pathologist ----

library(tidyverse)
df_compare <- readxl::read_xlsx("../Fig_codes/data/Fig4/GEM_Model_Percent_Correlation.xlsx")
lm <- lm(df_compare$Model_percent_0408 ~ df_compare$Human_percent, data = df_compare)
summary(lm)

ggplot(df_compare, aes(x = Human_percent, y = Model_percent_0408)) + 
  geom_point(size = 2.5, color = "#3c77af") + 
  geom_smooth(method = "lm", se = F, color = "black") + 
  labs(x = "Pathologist (%)", y = "AI (%)") + 
  theme_classic() + 
  theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
        plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
        legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='none',legend.text=element_text(size=14,hjust=0,face='bold'),
        axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
        axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
ggsave("figures/Fig4/0412_scatter_Patho_model_correlation.pdf", width = 3.3, height = 3)

scc.test <- cor.test(df_compare$Model_percent_0408, df_compare$Human_percent, method = "spearman")
scc.test

pcc.test <- cor.test(df_compare$Model_percent_0408, df_compare$Human_percent, method = "pearson")
pcc.test
