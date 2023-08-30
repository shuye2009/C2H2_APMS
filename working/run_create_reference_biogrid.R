# create gold standard PPI sets
setwd("C:/RSYNC/AP_MS_C2H2/C2H2_ZFP_data")
perlcmd <- "C:/RSYNC/Repositories/C2H2_APMS/working/create_references_biogrid.pl"

bait_files <- c("ZNF_baitfile_saint_C2H2_ZFP.tab",
               "ZNF_baitfile_saint_C2H2_ZFP_GFP.tab",
               "ZNF_baitfile_saint_C2H2_ZFP_TF_neg.tab")

prey_files <- c("ZNF_preyfile_saint_C2H2_ZFP.tab",
               "ZNF_preyfile_saint_C2H2_ZFP_GFP.tab",
               "ZNF_preyfile_saint_C2H2_ZFP_TF_neg.tab")
data_names <- c("all_neg", "GFP", "TF_neg")
bait_col <- 1 # perl use 0-based index
prey_col <- 0
id_col <- 0

for(i in seq_along(bait_files)){
  system(paste("perl", perlcmd, bait_files[i], prey_files[i], data_names[i], 
               bait_col, prey_col, id_col))
}

