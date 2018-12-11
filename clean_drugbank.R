# Required libraries: readxl,data.table,tidyverse

clean_drugbank <- function() {
  
  vocab <- read_csv("data/drugbank/drugbank_vocab.csv")
  vocab$substance <- paste0(vocab$`Common name`," | ",vocab$Synonyms)
  vocab <- vocab[,c("DrugBank ID","substance")]
  colnames(vocab) <- c("drugbank_id","substance")
  vocab$substance <- gsub(" \\| ",";",vocab$substance) 
  vocab <- separate_rows(vocab,substance,sep = ";")
  
  active <- read_csv("data/drugbank/drug_target_identifiers_all_pharmacologically_active_v5.1.1.csv")
  active <- active[,c("Gene Name","Drug IDs")]
  colnames(active) <- c("gene","drugbank_id")
  active <- separate_rows(active, drugbank_id)
  
  df <- merge(vocab,active,by = c("drugbank_id"),all.x = TRUE)
  
  df <- df[df$substance!="NA",]
  df$substance <- tolower(df$substance)
  df <- unique(df)
  
  return(df)
  
}