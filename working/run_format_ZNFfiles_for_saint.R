# format spectral counts file for SAINTexpress analysis

setwd("C:/RSYNC/AP_MS_C2H2/C2H2_ZFP_data")
perlcmd <- "C:/RSYNC/Repositories/C2H2_APMS/working/format_ZNFfiles_for_saint.pl"

protein_length_file <- "C:/RSYNC/AP_MS_C2H2/human_ids_length.txt";
spectral_count_file1 <- "C:/RSYNC/AP_MS_C2H2/C2H2_ZFP_data/1_All_C2H2_ZNF_Purifications_raw_data.txt"; #"20190118_AP_MS_data_for_SAINT_analysis.csv";#"160726_C2H2_ZNF_Purifications_fixed_plus_180826_New_purification.txt";
spectral_count_file2 <- "C:/RSYNC/AP_MS_C2H2/C2H2_ZFP_data/1_160726_TF_Negative_Control_Purifications.txt";
spectral_count_file3 <- "C:/RSYNC/AP_MS_C2H2/C2H2_ZFP_data/2_GFP_Purifications_raw_data.txt";

analysis <- "C2H2_ZFP";

remove_list <- "C:/RSYNC/AP_MS_C2H2/PPI_SAINT_2018/data180826/160816ListofproteinstoremovfromSAINTnetwork.txt";
update_list <- "C:/RSYNC/AP_MS_C2H2/human_name_update.txt";

system(paste("perl", perlcmd, 
             protein_length_file, 
             spectral_count_file1,
             spectral_count_file2, 
             spectral_count_file3, 
             analysis, 
             remove_list,
             update_list))