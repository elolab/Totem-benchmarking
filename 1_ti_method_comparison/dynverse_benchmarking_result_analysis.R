





data_set_info <- list()
for (file in list.files("~/linear_tree/"))
{
  dynobject <- readRDS(paste0("~/linear_tree/",file))
  data_set_info[[file]] <- c(dynobject$source,length(dynobject$milestone_ids),dynobject$trajectory_type)
  
}
data_set_info <- do.call(rbind,data_set_info)
rownames(data_set_info) <- gsub(".rds","",paste0("linear_tree_",rownames(data_set_info)))




library(ggplot2)

load("~/Totem_slingshot_results.RData")

ggplot(results,aes(x=type,y=overall,fill=type)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75)) +
  geom_boxplot(width=0.3)




filetable <- table(gsub("_single.rds|_type_1|_type_2|_type_3|_type_4|_type_5","",rownames(results)))

datasets_to_include <- names(filetable)[filetable==5]

rownamesresults <- rownames(results)[gsub("_single.rds|_type_1|_type_2|_type_3|_type_4|_type_5","",rownames(results)) %in% datasets_to_include]

results <- results[gsub("_single.rds|_type_1|_type_2|_type_3|_type_4|_type_5","",rownames(results)) %in% datasets_to_include,]

rownames(results) <- rownamesresults

ggplot(results,aes(x=type,y=correlation,fill=type)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75)) +
  geom_boxplot(width=0.3)

lapply(split(results$overall,results$type),summary)



file <- rownames(results)[results$type=="type_1_single"]
results_scshaper <- results[results$type=="type_1_single",]
rownames(results_scshaper) <- gsub("\\_type.*","",do.call(rbind,strsplit(file,"/"))[,3])
results_scshaper$file <- gsub("\\_type.*","",do.call(rbind,strsplit(file,"/"))[,3])
summary(results_scshaper$overall)
results_scshaper <- as.data.frame(results_scshaper)
results_scshaper$trajectory_type <- plyr::mapvalues(results_scshaper$file,from = rownames(data_set_info),to = data_set_info[,3])
results_scshaper$source <- plyr::mapvalues(results_scshaper$file,from = rownames(data_set_info),to = data_set_info[,1])


load("~/../Downloads/slingshotcustom_results.RData")
rownames(results) <- gsub(".rds","",do.call(rbind,strsplit(results$file,"/"))[,3])
files <- do.call(rbind,strsplit(datasets_to_include,"/"))[,3]
results <- results[do.call(rbind,strsplit(datasets_to_include,"/"))[,3],]
results$file <- files
results$overall <- (results$correlation)*(results$featureimp_wcor)*(results$F1_branches)*(results$him)
results$overall <- (results$overall)^(1/4)
summary(results$overall)
results_slingshot <- as.data.frame(results)
results_slingshot$trajectory_type <- plyr::mapvalues(results_slingshot$file,from = rownames(data_set_info),to = data_set_info[,3])
results_slingshot$source <- plyr::mapvalues(results_slingshot$file,from = rownames(data_set_info),to = data_set_info[,1])






load("~/../Downloads/Tinganew_seed_123456.RData")
rownames(results) <- gsub(".rds","",do.call(rbind,strsplit(results$file,"/"))[,3])
files <- do.call(rbind,strsplit(datasets_to_include,"/"))[,3]
results <- results[do.call(rbind,strsplit(datasets_to_include,"/"))[,3],]
results$file <- files
results$overall <- (results$correlation)*(results$featureimp_wcor)*(results$F1_branches)*(results$him)
results$overall <- (results$overall)^(1/4)
summary(results$overall)
results_Tinga <- as.data.frame(results)
results_Tinga$trajectory_type <- plyr::mapvalues(results_Tinga$file,from = rownames(data_set_info),to = data_set_info[,3])
results_Tinga$source <- plyr::mapvalues(results_Tinga$file,from = rownames(data_set_info),to = data_set_info[,1])



df <- results_slingshot
df$method <- "Slingshot"

df_ <- results_Tinga
df_$method <- "Tinga"
df <- rbind(df,df_)

df_ <- results_scshaper
df_$method <- "scShaper"
df_$type <- NULL
df <- rbind(df,df_)

df$synthetic <- grepl("synthetic",df$source)
df$linear <- df$trajectory_type=="linear"

df$method <- plyr::mapvalues(df$method,from = c("scShaper","Slingshot","Tinga"),
                             to = c("Totem","dynverse_Slingshot","Tinga"))

p1 <- ggplot(df,aes(x=method,y=correlation,fill=method)) +
  geom_violin(scale = "area") +
  geom_boxplot(width=0.1,fill="white") +
  ylab("Correlation") +
  xlab("Method") +
  theme_bw() +
  guides(fill=guide_legend(title="Method"))

p2 <- ggplot(df,aes(x=method,y=him,fill=method)) +
  geom_violin(scale = "area") +
  geom_boxplot(width=0.1,fill="white") +
  ylab("HIM") +
  xlab("Method") +
  theme_bw() +
  guides(fill=guide_legend(title="Method"))

p3 <- ggplot(df,aes(x=method,y=featureimp_wcor,fill=method)) +
  geom_violin(scale = "area") +
  geom_boxplot(width=0.1,fill="white") +
  ylab("Feature importance - weighted correlation") +
  xlab("Method") +
  theme_bw() +
  guides(fill=guide_legend(title="Method"))

p4 <- ggplot(df,aes(x=method,y=F1_branches,fill=method)) +
  geom_violin(scale = "area") +
  geom_boxplot(width=0.1,fill="white") +
  ylab("F1 branches") +
  xlab("Method") +
  theme_bw() +
  guides(fill=guide_legend(title="Method"))

p5 <- ggplot(df,aes(x=method,y=overall,fill=method)) +
  geom_violin(scale = "area") +
  geom_boxplot(width=0.1,fill="white") +
  ylab("Overall score") +
  xlab("Method") +
  theme_bw() +
  guides(fill=guide_legend(title="Method"))

cowplot::plot_grid(p1,p2,p3,p4,p5)



df_ <- reshape2::melt(df)
df_ <- df_[df_$variable %in% c("overall","correlation","him","F1_branches","featureimp_wcor"),]
df_$linear <- plyr::mapvalues(df_$linear,from = c("TRUE","FALSE"),c("linear","non-linear"))
df_$synthetic <- plyr::mapvalues(df_$synthetic,from = c("TRUE","FALSE"),c("synthetic","real"))

ggplot(df_,aes(x=method,y=value,fill=method)) +
  geom_violin(scale = "area") +
  geom_boxplot(width=0.1,fill="white") +
  ylab("Overall score") +
  xlab("Method") +
  theme_bw() +
  guides(fill=guide_legend(title="Method")) +
  facet_wrap(vars(variable))

p <- ggplot(df_,aes(x=method,y=value,fill=method)) +
  geom_violin(scale = "area") +
  geom_boxplot(width=0.1,fill="white") +
  ylab("Value") +
  xlab("Method") +
  theme_bw() +
  guides(fill=guide_legend(title="Method")) +
  facet_grid(linear~variable)

ggsave("~/totem_fig_1.pdf",plot = p,units = "in",width = 12,height = 5)



a1 <- df_[df_$variable=="overall" & df_$linear == "non-linear" & df_$method=="Totem",c("file","value")]
a2 <- df_[df_$variable=="overall" & df_$linear == "non-linear" & df_$method=="dynverse_Slingshot",c("file","value")]
a3 <- df_[df_$variable=="overall" & df_$linear == "non-linear" & df_$method=="Tinga",c("file","value")]
a1_ <- a1$value ; names(a1_) <- a1$file
a2_ <- a2$value ; names(a2_) <- a2$file
a3_ <- a3$value ; names(a3_) <- a3$file

wilcox.test(a1_,a2_[names(a1_)],paired = TRUE,alternative="two.sided")
wilcox.test(a1_,a3_[names(a1_)],paired = TRUE,alternative="two.sided")
# 
# wilcox.test(a1_,a2_[names(a1_)],paired = TRUE,alternative="greater")
# wilcox.test(a1_,a3_[names(a1_)],paired = TRUE,alternative="greater")
# 
# wilcox.test(a1_,a2_[names(a1_)],paired = TRUE,alternative="less")
# wilcox.test(a1_,a3_[names(a1_)],paired = TRUE,alternative="less")




a1 <- df_[df_$variable=="overall" & df_$linear == "linear" & df_$method=="Totem",c("file","value")]
a2 <- df_[df_$variable=="overall" & df_$linear == "linear" & df_$method=="dynverse_Slingshot",c("file","value")]
a3 <- df_[df_$variable=="overall" & df_$linear == "linear" & df_$method=="Tinga",c("file","value")]
a1_ <- a1$value ; names(a1_) <- a1$file
a2_ <- a2$value ; names(a2_) <- a2$file
a3_ <- a3$value ; names(a3_) <- a3$file

wilcox.test(a1_,a2_[names(a1_)],paired = TRUE,alternative="two.sided")
wilcox.test(a1_,a3_[names(a1_)],paired = TRUE,alternative="two.sided")

# wilcox.test(a1_,a2_[names(a1_)],paired = TRUE,alternative="greater")
# wilcox.test(a1_,a3_[names(a1_)],paired = TRUE,alternative="greater")
# 
# wilcox.test(a1_,a2_[names(a1_)],paired = TRUE,alternative="less")
# wilcox.test(a1_,a3_[names(a1_)],paired = TRUE,alternative="less")





p <- ggplot(df_,aes(x=method,y=value,fill=method)) +
  geom_violin(scale = "area") +
  geom_boxplot(width=0.1,fill="white") +
  ylab("Value") +
  xlab("Method") +
  theme_bw() +
  guides(fill=guide_legend(title="Method")) +
  facet_grid(synthetic~variable)

ggsave("~/totem_fig_1_2.pdf",plot = p,units = "in",width = 12,height = 5)


df_o <- df_[df_$variable=="overall",]

p <- ggplot(df_o,aes(x=method,y=value,fill=method)) +
  geom_violin(scale = "area") +
  geom_boxplot(width=0.1,fill="white") +
  ylab("Value") +
  xlab("Method") +
  theme_bw() +
  guides(fill=guide_legend(title="Method")) +
  facet_grid(trajectory_type~source)

ggsave("~/totem_fig_1_3.pdf",plot = p,units = "in",width = 12,height = 8)



ggplot(df,aes(x=method,y=overall,fill=method)) +
  geom_violin(scale = "area") +
  geom_boxplot(width=0.1,fill="white") +
  facet_wrap(vars(linear))



ggplot(df,aes(x=method,y=overall,fill=method)) +
  geom_violin(scale = "area") +
  geom_boxplot(width=0.1) +
  facet_wrap(vars(trajectory_type))


ggplot(df,aes(x=method,y=overall,fill=method)) +
  geom_violin(scale = "area") +
  geom_boxplot(width=0.1) +
  facet_wrap(vars(source))


ggplot(df,aes(x=method,y=overall,fill=method)) +
  geom_violin(scale = "area") +
  geom_boxplot(width=0.1) +
  facet_wrap(vars(synthetic))


ggplot(df,aes(x=method,y=overall,fill=method)) +
  geom_violin(scale = "area") +
  geom_boxplot(width=0.1) +
  facet_wrap(vars(linear))



ggplot(df,aes(x=method,y=overall,fill=method)) +
  geom_violin(scale = "area") +
  geom_boxplot(width=0.1,fill="white")

# facet_wrap(vars(trajectory_type))







ggplot(df,aes(x=method,y=overall,fill=method)) +
  geom_violin(scale = "area") +
  geom_boxplot(width=0.1) +
  facet_grid(trajectory_type~synthetic)
















