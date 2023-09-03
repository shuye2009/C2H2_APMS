wd <- "C:/RSYNC/ID_mapping"
perlcmd <- "C:/RSYNC/Repositories/C2H2_APMS/working/update_hgnc_gene_name_2019.pl"
setwd(wd)

infile <- "hgnc_gene_map_2023.txt"
outfile <- "hgnc_geneName_updates_2023.txt"

system(paste("perl", perlcmd, 
             infile, outfile))