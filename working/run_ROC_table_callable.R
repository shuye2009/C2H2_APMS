# perform ROC curve analysis and plot

library(dplyr)
library(ggplot2)
library(cowplot)

setwd("C:/RSYNC/AP_MS_C2H2/C2H2_ZFP_data")
perlcmd <- "C:/RSYNC/Repositories/C2H2_APMS/working/ROC_table_callable.pl"

score_files <- c("C2H2_ZFP_SAINT_score_list.txt",
                "C2H2_ZFP_SAINT_score_list_GFP.txt",
                "C2H2_ZFP_SAINT_score_list_TF_neg.txt") #comma seperated strings
score_columns <- rep(12, 3) #comma seperated integers, perl index is 0 based
sizes <- rep("small", 3) #comma seperated strings
gi_cols <- rep(0, 3) #single integer
baitfiles <- c("ZNF_baitfile_saint_C2H2_ZFP.tab",
              "ZNF_baitfile_saint_C2H2_ZFP_GFP.tab",
              "ZNF_baitfile_saint_C2H2_ZFP_TF_neg.tab") #single string
preyfiles <- c("ZNF_preyfile_saint_C2H2_ZFP.tab",
              "ZNF_preyfile_saint_C2H2_ZFP_GFP.tab",
              "ZNF_preyfile_saint_C2H2_ZFP_TF_neg.tab") #single string#single string
datasets <- c("all_neg", "GFP", "TF_neg") #single string

ROCs <- lapply(seq_along(score_files), function(i){
  system(paste("perl", perlcmd, score_files[i], score_columns[i], sizes[i],
               gi_cols[i], baitfiles[i], preyfiles[i], datasets[i]))
  res_file <- paste0(gsub("\\.txt", "", score_files[i]), "_GS_ROC_", sizes[i], ".tab")
  roc_table <- read.delim2(res_file, header = TRUE) %>%
    mutate(Background = datasets[i]) %>%
    mutate(SAINT_score = as.numeric(Cutoff),
           Sensitivity = as.numeric(sensitivity),
           Precision = as.numeric(precision),
           TPR = as.numeric(TPR),
           FPR = as.numeric(FPR))
})

rocs <- bind_rows(ROCs)

p1 <- ggplot(rocs, aes(x = FPR, y = TPR, color = Background)) +
  geom_point()
p1

p2 <- ggplot(rocs, aes(x = SAINT_score, y = Precision, color = Background)) +
  geom_point()
p2
p <- plot_grid(p1, p2, nrow = 2)
pdf("ROCs_plot.pdf", height = 10, width = 7)
print(p)
dev.off()
