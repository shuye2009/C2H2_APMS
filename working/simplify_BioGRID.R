## simplify BioGRID positive references
library(dplyr)
setwd("C:/RSYNC/ID_mapping")

mv_file <- "BIOGRID-MV-Physical-4.4.224.tab3.txt"
human_file <- "BIOGRID-ORGANISM-Homo_sapiens-4.4.224.tab3.txt"

mv <- read.delim2(mv_file, header = TRUE)
head(mv)
mv_simplified <- mv %>%
  filter(Organism.ID.Interactor.A == 9606 & Organism.ID.Interactor.B == 9606) %>%
  filter(Experimental.System.Type == "physical") %>%
  select(c(Official.Symbol.Interactor.A, Official.Symbol.Interactor.B, 
           SWISS.PROT.Accessions.Interactor.A, SWISS.PROT.Accessions.Interactor.B))
head(mv_simplified)
write.table(mv_simplified, "BIOGRID-MV-Physical-4.4.224.tab3_human.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

human <- read.delim2(human_file, header = TRUE)
head(human)
human_simplified <- human %>%
  filter(Organism.ID.Interactor.A == 9606 & Organism.ID.Interactor.B == 9606) %>%
  filter(Experimental.System.Type == "physical") %>%
  select(c(Official.Symbol.Interactor.A, Official.Symbol.Interactor.B, 
           SWISS.PROT.Accessions.Interactor.A, SWISS.PROT.Accessions.Interactor.B))
head(human_simplified)
write.table(human_simplified, "BIOGRID-ORGANISM-Homo_sapiens-4.4.224.tab3_human.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
