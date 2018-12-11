# Required libraries: readxl,data.table,tidyverse

clean_bnf <- function() {
  
  # Load BNF data ===============================================================
  
  df <- read.xlsx("data/Exposures.xlsx", sheetName = "IncludedExposures")
  
  # Restirct to antihypertensives ===============================================
  
  df <- df[df$Exposure=="Hypertension Drugs",]
  
  # Format dataframe ============================================================
  
  df <- df[,c("Sub.class","Drug.Substance.Name")]
  colnames(df) <- c("drug","substance")
  
  # Remove polypharmacy medines =================================================
  
  df <- df[!grepl("AND",df$drug,ignore.case = FALSE),]
  
  # Tidy drug substance information =============================================
  
  df$drug <- tolower(df$drug)
  df$substance <- tolower(df$substance)
  df <- df[!is.na(df$substance),]
  df <- df[!grepl("/",df$substance),]
  
  df[df$drug=="vasodilator antihypertensive drugs" & 
       df$substance=="sodium nitroprusside dihydrate",]$substance <- "nitroprusside"
  
  df$substance <- sub('(^\\w+)\\s.+','\\1',df$substance)
  
  
  df[df$drug=="thiazides and related diuretics" & 
       df$substance=="potassium",]$substance <- "potassium chloride"
  
  
  df <- unique(df)
  
  # Format drug names ===========================================================
  
  df$drug <- ifelse(df$drug=="beta-adrenoceptor blocking drugs",
                    "Beta-adrenoceptor blockers",df$drug)
  
  df$drug <- ifelse(df$drug=="renin inhibitors",
                    "Renin inhibitors",df$drug)
  
  df$drug <- ifelse(df$drug=="vasodilator antihypertensive drugs",
                    "Vasodilator antihypertensives",df$drug)
  
  df$drug <- ifelse(df$drug=="potassium-sparing diuretics and aldosterone antagonists",
                    "PSDs and aldosterone antagonists",df$drug)
  
  df$drug <- ifelse(df$drug=="calcium-channel blockers",
                    "Calcium channel blockers",df$drug)
  
  df$drug <- ifelse(df$drug=="centrally acting antihypertensive drugs",
                    "Centrally acting antihypertensives",df$drug)
  
  df$drug <- ifelse(df$drug=="thiazides and related diuretics",
                    "Thiazides and related diuretics",df$drug)
  
  df$drug <- ifelse(df$drug=="loop diuretics",
                    "Loop diuretics",df$drug)
  
  df$drug <- ifelse(df$drug=="angiotensin-ii receptor antagonists",
                    "Angiotensin-II receptor antagonists",df$drug)
  
  df$drug <- ifelse(df$drug=="angiotensin-converting enzyme inhibitors",
                    "Angiotensin converting enzyme inhibitors",df$drug)
  
  df$drug <- ifelse(df$drug=="alpha-adrenoceptor blocking drugs",
                    "Alpha-adrenoceptor blockers",df$drug)
  
  df$drug <- ifelse(df$drug=="adrenergic neurone blocking drugs",
                    "Adrenergic neurone blockers",df$drug)
  
  return(df)
  
}