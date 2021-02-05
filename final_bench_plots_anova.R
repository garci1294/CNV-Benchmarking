library(tidyverse)
library(ggpubr)
library(rstatix)
library(emmeans)

ICNV_A <- read.delim("~/Desktop/CNV_Benchmark/ICNV_AMP_STAT")
HB_A <- read.delim("~/Desktop/CNV_Benchmark/HB_AMP_STAT")
CS_A <- read.delim("~/Desktop/CNV_Benchmark/CS_AMP_STAT")
ICNV_D <- read.delim("~/Desktop/CNV_Benchmark/ICNV_DEL_STAT")
HB_D <- read.delim("~/Desktop/CNV_Benchmark/HB_DEL_STAT")
CS_D <- read.delim("~/Desktop/CNV_Benchmark/CS_DEL_STAT")

median_seconds <- c(ICNV_A$median, HB_A$median, CS_A$median, ICNV_D$median, HB_D$median, CS_D$median)
cnv_tool <- c(ICNV_A$tool, HB_A$tool, CS_A$tool, ICNV_D$tool, HB_D$tool, CS_D$tool)
ips <- c(ICNV_A$itr_sec, HB_A$itr_sec, CS_A$itr_sec, ICNV_D$itr_sec, HB_D$itr_sec, CS_D$itr_sec)
mem_allo <- c(ICNV_A$mem_alloc, HB_A$mem_alloc, CS_A$mem_alloc, ICNV_D$mem_alloc, HB_D$mem_alloc, CS_D$mem_alloc)
amp_del <- c(ICNV_A$supp, HB_A$supp, CS_A$supp, ICNV_D$supp, HB_D$supp, CS_D$supp)
total_time <- c(ICNV_A$total_time, HB_A$total_time, CS_A$total_time, ICNV_D$total_time, HB_D$total_time, CS_D$total_time)

median <- data.frame("Median"=median_seconds, "Type"=amp_del, "CNV_Tool"=cnv_tool)
itr_sec <- data.frame("Iteriation"=ips, "Type"=amp_del, "CNV_Tool"=cnv_tool)
memory <- data.frame("Memory"=mem_allo, "Type"=amp_del, "CNV_Tool"=cnv_tool)
total <- data.frame("Total"=total_time, "Type"=amp_del, "CNV_Tool"=cnv_tool)

median %>%
  group_by(Type, CNV_Tool) %>%
  get_summary_stats(Median, type = "mean_sd")

bxp_median <- ggboxplot(
  median, x = "Type", y = "Median", #legend = "",
  color = "CNV_Tool", palette = "jco",
  xlab = "", ylab = "Median (Seconds per Iteriation)",  legend.title = ""#, facet.by = "Type
)
bxp_median

median %>%
  group_by(Type, CNV_Tool) %>%
  identify_outliers(Median)

res.aov_median <- median %>% anova_test(Median ~ Type * CNV_Tool)
res.aov_median

model_median <- lm(Median ~ Type * CNV_Tool, data = median)
median %>%
  group_by(Type) %>%
  anova_test(Median ~ CNV_Tool, error = model_median)

# pairwise comparisons
pwc_median <- median %>% 
  group_by(Type) %>%
  emmeans_test(Median ~ CNV_Tool, p.adjust.method = "bonferroni") 
pwc_median

pwc_median <- pwc_median %>% add_xy_position(x = "Type")
bxp_median +
  stat_pvalue_manual(pwc_median) +
  labs(
    subtitle = get_test_label(res.aov_median, detailed = TRUE),
    caption = get_pwc_label(pwc_median)
  )
#-------------------------------------------------------------------------------

itr_sec %>%
  group_by(Type, CNV_Tool) %>%
  get_summary_stats(Iteriation, type = "mean_sd")

bxp_itr_sec <- ggboxplot(
  itr_sec, x = "Type", y = "Iteriation", #legend = "",
  color = "CNV_Tool", palette = "jco",
  xlab = "", ylab = "Iteriations per second",  legend.title = ""#, facet.by = "Type
)
bxp_itr_sec

itr_sec %>%
  group_by(Type, CNV_Tool) %>%
  identify_outliers(Iteriation)

res.aov_itr_sec <- itr_sec %>% anova_test(Iteriation ~ Type * CNV_Tool)
res.aov_itr_sec

model_itr_sec <- lm(Iteriation ~ Type * CNV_Tool, data = itr_sec)
itr_sec %>%
  group_by(Type) %>%
  anova_test(Iteriation ~ CNV_Tool, error = model_itr_sec)

# pairwise comparisons
pwc_itr_sec <- itr_sec %>% 
  group_by(Type) %>%
  emmeans_test(Iteriation ~ CNV_Tool, p.adjust.method = "bonferroni") 
pwc_itr_sec

pwc_itr_sec <- pwc_itr_sec %>% add_xy_position(x = "Type")
bxp_itr_sec +
  stat_pvalue_manual(pwc_itr_sec) +
  labs(
    subtitle = get_test_label(res.aov_itr_sec, detailed = TRUE),
    caption = get_pwc_label(pwc_itr_sec)
  )
#-------------------------------------------------------------------------------

memory %>%
  group_by(Type, CNV_Tool) %>%
  get_summary_stats(Memory, type = "mean_sd")

bxp_memory <- ggboxplot(
  memory, x = "Type", y = "Memory", #legend = "",
  color = "CNV_Tool", palette = "jco",
  xlab = "", ylab = "Allocated Memory/ 100 Iterations (Megabytes)",  legend.title = ""#, facet.by = "Type
)
bxp_memory

memory %>%
  group_by(Type, CNV_Tool) %>%
  identify_outliers(Memory)

res.aov_memory <- memory %>% anova_test(Memory ~ Type * CNV_Tool)
res.aov_memory

model_memory <- lm(Memory ~ Type * CNV_Tool, data = memory)
memory %>%
  group_by(Type) %>%
  anova_test(Memory ~ CNV_Tool, error = model_memory)

# pairwise comparisons
pwc_memory <- memory %>% 
  group_by(Type) %>%
  emmeans_test(Memory ~ CNV_Tool, p.adjust.method = "bonferroni") 
pwc_memory

pwc_memory <- pwc_memory %>% add_xy_position(x = "Type")
bxp_memory +
  stat_pvalue_manual(pwc_memory) +
  labs(
    subtitle = get_test_label(res.aov_memory, detailed = TRUE),
    caption = get_pwc_label(pwc_memory)
  )

#-------------------------------------------------------------------------------

total %>%
  group_by(Type, CNV_Tool) %>%
  get_summary_stats(Total, type = "mean_sd")

bxp_total <- ggboxplot(
  total, x = "Type", y = "Total", #legend = "",
  color = "CNV_Tool", palette = "jco",
  xlab = "", ylab = "Total time/ 100 Iteriations (minutes)", legend.title = ""#, facet.by = "Type
)
bxp_total

total %>%
  group_by(Type, CNV_Tool) %>%
  identify_outliers(Total)

res.aov_total <- total %>% anova_test(Total ~ Type * CNV_Tool)
res.aov_total

model_total <- lm(Total ~ Type * CNV_Tool, data = total)
total %>%
  group_by(Type) %>%
  anova_test(Total ~ CNV_Tool, error = model_total)

# pairwise comparisons
pwc_total <- total %>% 
  group_by(Type) %>%
  emmeans_test(Total ~ CNV_Tool, p.adjust.method = "bonferroni") 
pwc_total

pwc_total <- pwc_total %>% add_xy_position(x = "Type")
bxp_total +
  stat_pvalue_manual(pwc_total) +
  labs(
    subtitle = get_test_label(res.aov_total, detailed = TRUE),
    caption = get_pwc_label(pwc_total)
  )

