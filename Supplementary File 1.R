library(data.table);
library(dplyr);
library(ggplot2);
library(gridExtra);
library(openxlsx)
setwd("E:\\袤醱\\蜇翹")
###########Figure 6A~C#########################################################
boxplot<-read.xlsx("Supplementary Table 1.xlsx",colNames = T)
boxplot$subgroup <- as.factor(boxplot$subgroup)
median <- data.frame(subgroup=seq(1, 10, 1), median=tapply(boxplot$CH12_log_H3K79me2_average, boxplot$subgroup, median))
variance <- data.frame(subgroup=seq(1, 10, 1), variance=tapply(boxplot$log2FoldChange, boxplot$subgroup, var))
boxplot <- merge(boxplot, median, all=T, by="subgroup")
#######Figure 1D###############
x <- factor(boxplot$subgroup)
y<- boxplot$CH12_log_H3K79me2_average
df <- boxplot
fill <- boxplot$median

ggplot(df, aes(x=x, y=y, fill=fill)) + 
  stat_boxplot(geom="errorbar", width=0.2, lwd=0.4) +
  geom_boxplot(width=0.6, position=position_dodge(0.8), show.legend=F) +
  scale_fill_gradient(low = "#00AFC1", high = "#006778") +
  xlab("Subgroup") +
  ylab("CH12 log2(H3K79me2)") +
  theme_bw() + 
  theme(text = element_text(family = "serif"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14),
        axis.ticks.length=unit(-0.2, "cm"),
        panel.background=element_rect(fill='transparent'),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.9)
        )
ggsave("Figure_6A.pdf", width = 7, height = 4)

x_variance <- factor(variance$subgroup)
y_variance <- variance$variance
df <- variance
fill <- variance$variance

ggplot(df, aes(x=x_variance, y=y_variance, fill=fill)) +
  geom_bar(stat="identity", width=0.8, show.legend=F) +
  xlab("Subgroup") +
  ylab("Variance") +
  scale_y_continuous(expand=c(0.008,0), limits= c(0, 35)) +
  scale_fill_gradient(low="#548CA8", high="#334257") +
  theme_bw() + 
  theme(text = element_text(family = "serif"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14),
        axis.ticks.length.y=unit(-0.2, "cm"),
        axis.ticks.length.x=unit(-0.05, "cm"),
        panel.background=element_rect(fill='transparent'),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.9)
        )
ggsave("Figure_6B.pdf", width = 7, height = 4)
  
df <- merge(median,variance,by='subgroup')
x <- df$median
y <- df$variance
老 <- cor.test(x,y,method = c("spearman"),exact=FALSE)$estimate
老
pvalue<-cor.test(x,y,method = c("spearman"),exact=FALSE)$p.value
pvalue
ggplot(df, aes(x=x, y=y))+
  geom_point(size=3, shape=18) +
  geom_text(aes(x=0, y=33),label=expression(paste(rho, " = -0.93", ", p =7.32e-06")),size=7, family="serif", face="bold")+
  geom_smooth(method=lm, color="red",se=F, size=1.3) +
  xlab("CH12 log2(H3K79me2) median") +
  ylab("Variance") +
  scale_x_continuous(breaks=seq(-7, 2.5, 2), limits=c(-7, 2.5), expand=c(0.005, 0)) +
  scale_y_continuous(breaks = seq(5, 35, 10), limits = c(5, 35)) +
  theme_bw() + 
  theme(text = element_text(family = "serif"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14),
        axis.ticks.length=unit(-0.2, "cm"),
        panel.background=element_rect(fill='transparent'),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.9)
        )
ggsave("Figure_6C.pdf", width = 7, height = 4)

################Figure 6D~F###########################################
boxplot<-read.xlsx("Supplementary Table 2.xlsx",colNames = T)
boxplot$subgroup <- as.factor(boxplot$subgroup)
median <- data.frame(subgroup=seq(1, 10, 1), median=tapply(boxplot$WT_log_H3K79me2, boxplot$subgroup, median))
variance <- data.frame(subgroup=seq(1, 10, 1), variance=tapply(boxplot$log2FoldChange, boxplot$subgroup, var))
boxplot <- merge(boxplot, median, all=T, by="subgroup")

x_median <- factor(boxplot$subgroup)
y_median <- boxplot$WT_log_H3K79me2
df <- boxplot
fill <- boxplot$median

ggplot(df, aes(x=x_median, y=y_median, fill=fill)) + 
  stat_boxplot(geom="errorbar", width=0.2, lwd=0.4) +
  geom_boxplot(width=0.6, position=position_dodge(0.8), show.legend=F) +
  scale_fill_gradient(low = "#00AFC1", high = "#006778") +
  xlab("Subgroup") +
  ylab("WT log2(H3K79me2)") +
  theme_bw() + 
  theme(text = element_text(family = "serif"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14),
        axis.ticks.length=unit(-0.2, "cm"),
        panel.background=element_rect(fill='transparent'),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.9)
        )
ggsave("Figure_6D.pdf", width = 7, height = 4)

x_variance <- factor(variance$subgroup)
y_variance <- variance$variance
df <- variance
fill <- variance$variance

ggplot(df, aes(x=x_variance, y=y_variance, fill=fill)) +
  geom_bar(stat="identity", width=0.8, show.legend=F) +
  xlab("Subgroup") +
  ylab("Variance") +
  scale_y_continuous(expand=c(0.008,0), limits= c(0, 35)) +
  scale_fill_gradient(low="#548CA8", high="#334257") +
  theme_bw() + 
  theme(text = element_text(family = "serif"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14),
        axis.ticks.length.y=unit(-0.2, "cm"),
        axis.ticks.length.x=unit(-0.05, "cm"),
        panel.background=element_rect(fill='transparent'),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.9)
        )
ggsave("Figure_6E.pdf", width = 7, height = 4)
  

df <- merge(median,variance,by='subgroup')
x <- df$median
y <- df$variance
老 <- cor.test(x,y,method = c("spearman"),exact=FALSE)$estimate
老
pvalue<-cor.test(x,y,method = c("spearman"),exact=FALSE)$p.value
pvalue
ggplot(df, aes(x=x, y=y))+
  geom_text(aes(x=0.5, y=33),label=expression(paste(rho, " = -1", ", p =1.70e-61")),size=7, family="serif", face="bold")+
  geom_point(size=3, shape=18) +
  geom_smooth(method=lm, color="red",se=F, size=1.3) +
  xlab("WT log2(H3K79me2) median") +
  ylab("Variance") +
  scale_x_continuous(breaks=seq(-4, 2, 2), limits=c(-4, 2), expand=c(0.005, 0)) +
  scale_y_continuous(breaks = seq(5, 35, 10), limits = c(5, 35)) +
  theme_bw() + 
  theme(text = element_text(family = "serif"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14),
        axis.ticks.length=unit(-0.2, "cm"),
        panel.background=element_rect(fill='transparent'),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.9)
        )
ggsave("Figure_6F.pdf", width = 7, height = 4)

###########Figure 6G~I#########################################################
boxplot<-read.xlsx("Supplementary Table 3.xlsx",colNames = T)
boxplot$subgroup <- as.factor(boxplot$subgroup)
median <- data.frame(subgroup=seq(1, 10, 1), median=tapply(boxplot$IMR90_log_H3K79me2_average, boxplot$subgroup, median))
variance <- data.frame(subgroup=seq(1, 10, 1), variance=tapply(boxplot$log2FoldChange, boxplot$subgroup, var))
boxplot <- merge(boxplot, median, all=T, by="subgroup")
x <- factor(boxplot$subgroup)
y<- boxplot$IMR90_log_H3K79me2_average
df <- boxplot
fill <- boxplot$median

ggplot(df, aes(x=x, y=y, fill=fill)) + 
  stat_boxplot(geom="errorbar", width=0.2, lwd=0.4) +
  geom_boxplot(width=0.6, position=position_dodge(0.8), show.legend=F) +
  scale_fill_gradient(low = "#00AFC1", high = "#006778") +
  xlab("Subgroup") +
  ylab("IMR90 log2(H3K79me2)") +
  theme_bw() + 
  theme(text = element_text(family = "serif"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14),
        axis.ticks.length=unit(-0.2, "cm"),
        panel.background=element_rect(fill='transparent'),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.9)
        )
ggsave("Figure_6G.pdf", width = 7, height = 4)

x_variance <- factor(variance$subgroup)
y_variance <- variance$variance
df <- variance
fill <- variance$variance

ggplot(df, aes(x=x_variance, y=y_variance, fill=fill)) +
  geom_bar(stat="identity", width=0.8, show.legend=F) +
  xlab("Subgroup") +
  ylab("Variance") +
  scale_y_continuous(expand=c(0.008,0), limits= c(0, 40)) +
  scale_fill_gradient(low="#548CA8", high="#334257") +
  theme_bw() + 
  theme(text = element_text(family = "serif"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14),
        axis.ticks.length.y=unit(-0.2, "cm"),
        axis.ticks.length.x=unit(-0.05, "cm"),
        panel.background=element_rect(fill='transparent'),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.9)
        )
ggsave("Figure_6H.pdf", width = 7, height = 4)
  

df <- merge(median,variance,by='subgroup')
x <- df$median
y <- df$variance
老 <- cor.test(x,y,method = c("spearman"),exact=FALSE)$estimate
老
pvalue<-cor.test(x,y,method = c("spearman"),exact=FALSE)$p.value
pvalue
ggplot(df, aes(x=x, y=y))+
  geom_point(size=3, shape=18) +
  geom_text(aes(x=0, y=33),label=expression(paste(rho, " = -0.82", ", p =3.81e-03")),size=7, family="serif", face="bold")+
  geom_smooth(method=lm, color="red",se=F, size=1.3) +
  xlab("IMR90 log2(H3K79me2) median") +
  ylab("Variance") +
  scale_x_continuous(breaks=seq(-5, 1.5, 1), limits=c(-5, 1.5), expand=c(0.005, 0)) +
  scale_y_continuous(breaks = seq(10, 38, 10), limits = c(10, 38)) +
  theme_bw() + 
  theme(text = element_text(family = "serif"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14),
        axis.ticks.length=unit(-0.2, "cm"),
        panel.background=element_rect(fill='transparent'),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.9)
        )
ggsave("Figure_6I.pdf", width = 7, height = 4)




################Figure 6J~L###########################################
boxplot<-read.xlsx("Supplementary Table 4.xlsx",colNames = T)
boxplot$subgroup <- as.factor(boxplot$subgroup)
median <- data.frame(subgroup=seq(1, 10, 1), median=tapply(boxplot$hepatocyte_log_H3K79me2, boxplot$subgroup, median))
variance <- data.frame(subgroup=seq(1, 10, 1), variance=tapply(boxplot$log2FoldChange, boxplot$subgroup, var))
boxplot <- merge(boxplot, median, all=T, by="subgroup")
x_median <- factor(boxplot$subgroup)
y_median <- boxplot$hepatocyte_log_H3K79me2
df <- boxplot
fill <- boxplot$median

ggplot(df, aes(x=x_median, y=y_median, fill=fill)) + 
  stat_boxplot(geom="errorbar", width=0.2, lwd=0.4) +
  geom_boxplot(width=0.6, position=position_dodge(0.8), show.legend=F) +
  scale_fill_gradient(low = "#00AFC1", high = "#006778") +
  xlab("Subgroup") +
  ylab("hepatocyte log2(H3K79me2)") +
  theme_bw() + 
  theme(text = element_text(family = "serif"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14),
        axis.ticks.length=unit(-0.2, "cm"),
        panel.background=element_rect(fill='transparent'),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.9)
        )
ggsave("Figure_6J.pdf", width = 7, height = 4)

x_variance <- factor(variance$subgroup)
y_variance <- variance$variance
df <- variance
fill <- variance$variance

ggplot(df, aes(x=x_variance, y=y_variance, fill=fill)) +
  geom_bar(stat="identity", width=0.8, show.legend=F) +
  xlab("Subgroup") +
  ylab("Variance") +
  scale_y_continuous(expand=c(0.008,0), limits= c(0, 35)) +
  scale_fill_gradient(low="#548CA8", high="#334257") +
  theme_bw() + 
  theme(text = element_text(family = "serif"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14),
        axis.ticks.length.y=unit(-0.2, "cm"),
        axis.ticks.length.x=unit(-0.05, "cm"),
        panel.background=element_rect(fill='transparent'),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.9)
        )
ggsave("Figure_6K.pdf", width = 7, height = 4)
  

df <- merge(median,variance,by='subgroup')
x <- df$median
y <- df$variance
老 <- cor.test(x,y,method = c("spearman"),exact=FALSE)$estimate
老
pvalue<-cor.test(x,y,method = c("spearman"),exact=FALSE)$p.value
pvalue
ggplot(df, aes(x=x, y=y))+
  geom_point(size=3, shape=18) +
  geom_text(aes(x=-0.6, y=30),label=expression(paste(rho, " = -0.98", ", p =1.47e-6")),size=7, family="serif", face="bold")+
  geom_smooth(method=lm, color="red",se=F, size=1.3) +
  xlab("IMR90 log2(H3K79me2) median") +
  ylab("Variance") +
  scale_x_continuous(breaks=seq(-3, 0.5, 1), limits=c(-3, 0.5), expand=c(0.005, 0)) +
  scale_y_continuous(breaks = seq(10, 32, 10), limits = c(10, 32)) +
  theme_bw() + 
  theme(text = element_text(family = "serif"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14),
        axis.ticks.length=unit(-0.2, "cm"),
        panel.background=element_rect(fill='transparent'),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.9)
        )
ggsave("Figure_6L.pdf", width = 7, height = 4)
