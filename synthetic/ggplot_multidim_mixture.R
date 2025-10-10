library(latex2exp)
library(tidyverse)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(scales)

# helper: load RData, normalize by ls_bs, return long df with dimension label
process_res2 <- function(file, d){
  load(file)  # loads res2
  tmp <- res2 / res2[,6]   # ratio wrt ls_bs
  colnames(tmp)[1] <- "alpha"
  tmp <- data.frame(tmp)
  tmp[,1] <- res2[,1]
  df <- tmp %>%
    gather(key = "Type", value = "value", -alpha)
  df[df=="forster_poly"] <- "FW_poly"
  df[df=="forster_bs"]   <- "FW_bs"
  df[df=="forster_ns"]   <- "FW_ns"
  df$dimension <- paste0("d = ", d)
  return(df)
}

# function to build plot for a given n across all d
build_plot <- function(n){
  dims <- c(1,5,10,50)
  dfs <- lapply(dims, function(d){
    fname <- paste0("mixture_forster_ls_", n, "_d", d, ".RData")
    process_res2(fname, d)
  })
  df_all <- bind_rows(dfs)
  
  ggplot(df_all, aes(x = alpha, y = value, color = Type, linetype = Type, shape = Type)) +
    geom_line(linewidth=2) +
    geom_point(size=5) +
    facet_wrap(~dimension, ncol=2) +
    scale_linetype_manual(values = c(rep("solid",3), rep("dashed",3))) +
    scale_color_manual(values = c(rep(hue_pal()(6)[1:3],2))) +
    ylab("Ratio of MSE") +
    xlab(TeX(paste0("($n=", n, "$)"))) +
    theme(strip.text = element_text(size=30, face="bold"),
          axis.ticks.length=unit(.4,"cm"),
          legend.position="top",
          legend.title = element_text(size=30, face="bold"),
          legend.text  = element_text(size=25),
          strip.background = element_blank(),
          strip.placement = "outside",
          axis.title.x = element_text(size=30, face="bold"),
          axis.text.x  = element_text(face="bold", size=25),
          axis.text.y  = element_text(face="bold", size=25),
          axis.title.y = element_text(face="bold", size=30))
}

# build plots for n=2000 and n=400
final1 <- build_plot(2000)


dims <- c(5,10,50)
dfs <- lapply(dims, function(d){
  fname <- paste0("mixture_forster_ls_", 2000, "_d", d, ".RData")
  process_res2(fname, d)
})
df_all <- bind_rows(dfs)

# summarize by taking mean over alpha values
summary_table <- df_all %>%
  group_by(dimension, Type) %>%
  summarise(mean_ratio = mean(value, na.rm=TRUE), .groups="drop") %>%
  spread(Type, mean_ratio)

print(summary_table)

