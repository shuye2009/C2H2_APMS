## Run SAINT command
setwd("C:/RSYNC/AP_MS_C2H2/C2H2_ZFP_data")

run1 <- system("SAINTexpress-spc \
       ZNF_interactionfile_saint_C2H2_ZFP.tab \
       ZNF_preyfile_saint_C2H2_ZFP.tab \
       ZNF_baitfile_saint_C2H2_ZFP.tab")
system("mv list.txt C2H2_ZFP_SAINT_score_list.txt")
run2 <- system("SAINTexpress-spc \
       ZNF_interactionfile_saint_C2H2_ZFP_GFP.tab \
       ZNF_preyfile_saint_C2H2_ZFP_GFP.tab \
       ZNF_baitfile_saint_C2H2_ZFP_GFP.tab")
system("mv list.txt C2H2_ZFP_SAINT_score_list_GFP.txt")
run3 <- system("SAINTexpress-spc \
       ZNF_interactionfile_saint_C2H2_ZFP_TF_neg.tab \
       ZNF_preyfile_saint_C2H2_ZFP_TF_neg.tab \
       ZNF_baitfile_saint_C2H2_ZFP_TF_neg.tab")
system("mv list.txt C2H2_ZFP_SAINT_score_list_TF_neg.txt")
