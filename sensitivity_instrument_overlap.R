# LOAD DATA ===================================================================

df <- data.frame(read_csv("data/mr_sbp_ad.csv"))
df <- df[df$level %in% c("Drug","Overall") & df$methodtype=="main",c("leveli","snps.mr")]
df <- df %>% 
  mutate(snp = strsplit(as.character(snps.mr), ";")) %>% 
  unnest(snp)

larsson <- read.xlsx("data/larsson.xlsx",sheetName = "Sheet1",
                     stringsAsFactors=FALSE)

larsson$snp <- gsub(" ","",larsson$snp)

ostergaard <- read.xlsx("data/ostergaard.xlsx",sheetName = "Sheet1",
                        stringsAsFactors=FALSE)

gill <- read.xlsx("data/gill.xlsx",sheetName = "Sheet1",
                 stringsAsFactors=FALSE)

# FIND SNP PROXIES IN OTHER STUDIES ===========================================

df$snps.mr <- NULL
df$snp.match <- NA
df$source.match <- NA
df$rsq <- NA
df$error <- FALSE

other_studies <- data.frame(snp = c(larsson$snp,ostergaard$snp,gill$SNP),
                 source = c(rep("Larsson",length(larsson$snp)),
                            rep("Ostergaard",length(ostergaard$snp)),
                            rep("Gill",length(gill$SNP))),
                 stringsAsFactors = FALSE)

for (i in c(1:nrow(df))) { 
  
  d <- tryCatch(get_proxies(query = df$snp[i]), error=function(e) NULL)
  d <- unique(d[,c("ID","R.squared")])
  
  if (is.null(d)) {
    df$error[i] <- TRUE
  } else {
    df$snp.match[i] <- paste0(d[d$ID %in% other_studies$snp,]$ID, collapse = ";")
    df$source.match[i] <- paste0(other_studies[other_studies$snp %in% (d[d$ID %in% other_studies$snp,]$ID),2], collapse = ";")
    df$rsq[i] <- paste0(sprintf("%.3f",(d[d$ID %in% other_studies$snp,]$R.squared)), collapse = ";")
  }
  
}

write.csv(df,file="output/instrument_overlap.csv",row.names = FALSE,na="")

df$error <- NULL
df <- df[!is.na(df$snp.match) & df$snp.match!="",]

df <- df %>% 
  mutate(snp.match = strsplit(as.character(snp.match), ";"),
         source.match = strsplit(as.character(source.match), ";"),
         rsq = strsplit(as.character(rsq), ";")) %>% 
  unnest(snp.match,source.match,rsq)

write.xlsx(data.frame(df), file="output/Supplementary_Tables_AntihypertensivesMR.xlsx", 
           sheetName="ST6",append=TRUE, row.names = FALSE, showNA = FALSE)
