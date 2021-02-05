library(ROCit)

CS_AMP_ROC <- read.delim("~/Desktop/CNV_Benchmark/CNV_ROC/CNV_AMP_ROC.txt")
CS_DEL_ROC <- read.delim("~/Desktop/CNV_Benchmark/CNV_ROC/CNV_DEL_ROC.txt")

# ================================================================================================
roc_empirical_AMP_CS <- rocit(score = CS_AMP_ROC$CS_amp_mb, class = CS_AMP_ROC$CS_amp_dtest,
                       negref = "-")

ciAUC(roc_empirical_AMP_CS)

plot(roc_empirical_AMP_CS, col = c(2,4), YIndex = T, values = F)


legend("topleft", "CONICSmat Amplification ROC")

legend("bottomright", c(    "Empirical ROC (CL = 95%)",
                       "-----------------------------------",
                        "estimated AUC : 0.492767123287671",
                        "lower = 0.409329055163602", 
                        "upper = 0.576205191411741"))

# ================================================================================================
roc_empirical_AMP_HB <- rocit(score = CS_AMP_ROC$HB_amp_mb, class = CS_AMP_ROC$HB_amp_dtest,
                            negref = "-")

ciAUC(roc_empirical_AMP_HB)

plot(roc_empirical_AMP_HB, col = c(2,4), YIndex = T, values = F)

legend("topleft", "HoneyBadger Amplification ROC")

legend("bottomright", c(     "Empirical ROC (CL = 95%)",
                       "-----------------------------------",
                        "estimated AUC : 0.706937095825985",
                        "lower = 0.55876886941119", 
                        "upper = 0.85510532224078"))

# ================================================================================================
roc_empirical_AMP_ICNV <- rocit(score = CS_AMP_ROC$ICNV_amp_mb, class = CS_AMP_ROC$INCV_amp_dtest,
                       negref = "-")

ciAUC(roc_empirical_AMP_ICNV)

plot(roc_empirical_AMP_ICNV, col = c(2,4), YIndex = T, values = F)

legend("topleft", "InferCNV Amplification ROC")

legend("bottomright", c(     "Empirical ROC (CL = 95%)",
                      "-----------------------------------",
                        "estimated AUC : 0.733484946935239",
                        "lower = 0.645334666047747", 
                        "upper = 0.821635227822732"))

# ================================================================================================
roc_empirical_DEL_CS <- rocit(score = CS_DEL_ROC$CS_del_mb, class = CS_DEL_ROC$CS_del_dtest,
                              negref = "-")

ciAUC(roc_empirical_DEL_CS)

plot(roc_empirical_DEL_CS, col = c(2,4), YIndex = T, values = F)

legend("topleft", "CONICSmat Deletion ROC")

legend("bottomright", c(    "Empirical ROC (CL = 95%)",
                            "-----------------------------------",
                            "estimated AUC : 0.628127696289905",
                            "lower = 0.546999080996858", 
                            "upper = 0.709256311582952"))

# ================================================================================================
roc_empirical_DEL_HB <- rocit(score = CS_DEL_ROC$HB_del_mb, class = CS_DEL_ROC$HB_del_dtest,
                              negref = "-")

ciAUC(roc_empirical_DEL_HB)

plot(roc_empirical_DEL_HB, col = c(2,4), YIndex = T, values = F)

legend("topleft", "HoneyBadger Deletion ROC")

legend("bottomright", c(     "Empirical ROC (CL = 95%)",
                             "-----------------------------------",
                             "estimated AUC : 0.758098006644518",
                             "lower = 0.692256209724678", 
                             "upper = 0.823939803564358"))

# ================================================================================================
roc_empirical_DEL_ICNV <- rocit(score = CS_DEL_ROC$ICNV_del_mb, class = CS_DEL_ROC$INCV_del_dtest,
                                negref = "-")

ciAUC(roc_empirical_DEL_ICNV)

plot(roc_empirical_DEL_ICNV, col = c(2,4), YIndex = T, values = F)

legend("topleft", "InferCNV Deletion ROC")

legend("bottomright", c(     "Empirical ROC (CL = 95%)",
                             "-----------------------------------",
                             "estimated AUC : 0.744393687707641",
                             "lower = 0.673800231013275", 
                             "upper = 0.814987144402007"))

# ================================================================================================


