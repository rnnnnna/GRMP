library("dplyr")
library("ggpubr")
library("tidyverse")
library(agricolae)
library(car)
library(reshape2)
library(dplyr)

#Source data : data_fig .1a.txt.
yield = read.table("data_fig .1a.txt.", header=T, sep="\t")

ggqqplot(yield$Yield)
ggdensity(yield$Yield,
          xlab = "yield")
shapiro.test(yield$Yield)#p>0.05
leveneTest(Yield~ group, data = yield)#p>0.05
levels(yield$group)
group_by(yield, group) %>%
  summarise(
    count = n(),
    mean = mean(Yield, na.rm = TRUE),
    sd = sd(Yield, na.rm = TRUE)
  )

res.aov <- aov(Yield ~ group, data = yield)
summary(res.aov)#p<0.05  
Tukey_HSD = TukeyHSD(res.aov, ordered = TRUE, conf.level = 0.95)
Tukey_HSD_table = as.data.frame(Tukey_HSD$group)
Tukey_HSD_table
write.table (Tukey_HSD_table, file ="yield_Tukey_HSD_table.txt", row.names = TRUE,col.names = TRUE, sep="\t")

out <- LSD.test(res.aov,"group", p.adj="none")
out
out$group
stat = out$groups
yield$stat=stat[as.character(yield$group),]$groups
yield$stat
yield[,"Yield"]
max=max(yield[,c("Yield")])
min=min(yield[,c("Yield")])
x = yield[,c("group","Yield")]
y = x %>% group_by(group) %>% summarise_(Max=paste('max(',"Yield",')',sep=""))
y=as.data.frame(y)
rownames(y)=y$group
yield$y=y[as.character(yield$group),]$Max + (max-min)*0.05
yield
theme_set(theme_classic())
p = ggplot(yield, aes(x=group, y=yield[["Yield"]], color=group)) +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
  labs(x="Groups", y=paste("Yield(kg/ha)")) + theme_classic() +
  geom_text(data=yield, aes(x=group, y=y, color=group, label= stat)) +
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)
if (length(unique(yield$group))>3){
  p=p+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))
}
color = c("#00BFC4", "#F9766E", "#00BFC4", "#F9766E", "#00BFC4", "#F9766E" )
p = p + scale_color_manual(values = color)
ggsave(paste(output_dir,"Fig .1a", ".pdf", sep=""), p, width = 5, height = 3)


#Source data : data_fig .1b.txt.
diameter = read.table("data_fig .1b.txt", header=T, sep="\t")
str(diameter)
ggqqplot(diameter$BulbDiameter)
ggdensity(diameter$BulbDiameter, xlab = "diameter")
shapiro.test(diameter$BulbDiameter)#p<0.05
leveneTest(BulbDiameter ~ group, data = diameter)#p<0.05
levels(diameter$group)
library(dplyr)
group_by(diameter, group) %>%
  summarise(
    count = n(),
    mean = mean(BulbDiameter, na.rm = TRUE),
    sd = sd(BulbDiameter, na.rm = TRUE)
  )

res.kru <- kruskal.test(BulbDiameter ~ group, data = diameter)
res.kru  #p<0.05

library(DescTools)
DunnTest(BulbDiameter~group, diameter)
out<-kruskal(diameter$BulbDiameter,diameter$group,group=TRUE,p.adj="bonferroni")

out$group
stat = out$groups
diameter$stat=stat[as.character(diameter$group),]$groups
diameter$stat
diameter[,"BulbDiameter"]
max=max(diameter[,c("BulbDiameter")])
min=min(diameter[,c("BulbDiameter")])
x = diameter[,c("group","BulbDiameter")]
y = x %>% group_by(group) %>% summarise_(Max=paste('max(',"BulbDiameter",')',sep=""))
y=as.data.frame(y)
rownames(y)=y$group
diameter$y=y[as.character(diameter$group),]$Max + (max-min)*0.05
diameter
theme_set(theme_classic())
p = ggplot(diameter, aes(x=group, y=diameter[["BulbDiameter"]], color=group)) +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
  labs(x="Groups", y=paste("Bulb Diameter(cm)")) + theme_classic() +
  geom_text(data=diameter, aes(x=group, y=y, color=group, label= stat)) +
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)
if (length(unique(diameter$group))>3){
  p=p+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))
}
color = c("#00BFC4", "#F9766E", "#00BFC4", "#F9766E", "#00BFC4", "#F9766E" )
p = p + scale_color_manual(values = color)
ggsave(paste(output_dir,"Fig 1.b", ".pdf", sep=""), p, width = 5, height = 3)


#Source data : data_fig .7.txt.
len = read.table("data_fig .7", header=T, sep="\t")
len = len %>% gather("group","length") %>%drop_na(.) 
table(len$group)

ggqqplot(len$length)
ggdensity(len$length, xlab = "len")
shapiro.test(len$length)#p<0.05
leveneTest(length ~ group, data = len)#p<0.05
len$group=factor(len$group,levels = c('CK','B','M','M+B'))
levels(len$group)

group_by(len, group) %>%
  summarise(
    count = n(),
    mean = mean(length, na.rm = TRUE),
    sd = sd(length, na.rm = TRUE)
  )
res.kru <- kruskal.test(length ~ group, data = len)
res.kru  #p<0.05
library(DescTools)
DunnTest(length~group, len)

out<-kruskal(len$length,len$group,group=TRUE,p.adj="bonferroni")
out$group
stat = out$groups
len$stat=stat[as.character(len$group),]$groups
len$stat
len[,"length"]
max=max(len[,c("length")])
min=min(len[,c("length")])
x = len[,c("group","length")]
y = x %>% group_by(group) %>% summarise_(Max=paste('max(',"length",')',sep=""))
y=as.data.frame(y)
rownames(y)=y$group
len$y=y[as.character(len$group),]$Max + (max-min)*0.05
theme_set(theme_classic())
p = ggplot(len, aes(x=group, y=len[["length"]], color=group)) +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
  labs(x="Groups", y=paste("Length(cm)")) + theme_classic() +
  geom_text(data=len, aes(x=group, y=y, color=group, label= stat)) +
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)

if (length(unique(len$group))>3){
  p=p+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))
}
color = c(  "green3", "#5F7FC7", "#D14285","gold","#DA5724", "orange","#508578", "#CD9BCD", "red3",
            "#AD6F3B", "#673770", "#652926", "#CBD588",
            "#8569D5" )
p = p + scale_color_manual(values = color)
ggsave(paste(output_dir,"Fig. 7", ".pdf", sep=""), p, width = 5, height = 3)
