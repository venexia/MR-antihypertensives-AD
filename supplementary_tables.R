# Supplementary Table 2 =======================================================

input_exp <- fread("data/input_exp.csv",
                   stringsAsFactors = FALSE,
                   data.table = FALSE)

input_exp <- input_exp[,c("gene","tissue","keep")]

colnames(input_exp) <- c("gene","tissue","ind.ins")

mr_exp_sbp <- fread("data/mr_exp_sbp.csv",
                    stringsAsFactors = FALSE,
                    data.table = FALSE)

st2 <- merge(input_exp,mr_exp_sbp,all.x = TRUE, by = c("gene","tissue"))

st2$ind.mr <- ifelse(!is.na(st2$beta.mr),TRUE,FALSE)

st2$ind.sbpev <- ifelse(!is.na(st2$beta.mr) & sign(st2$lci.mr)==sign(st2$uci.mr),TRUE,FALSE)

st2 <- st2[,c(1:2,11,19:20,24:25,21:23,6,9:10,3,26:27)]

colnames(st2) <- c(colnames(st2)[1:2],"SNP",colnames(st2)[4:7],
                   "beta.expression","se.expression","pval.expression",
                   colnames(st2)[11:16])
  
write.xlsx(data.frame(st2), file="output/Supplementary_Tables_AntihypertensivesMR.xlsx", 
           sheetName="ST2",append=TRUE, row.names = FALSE, showNA = FALSE)

# Supplementary Table 3 ======================================================= 

st3 <- fread("data/input_sbp.csv",
                   stringsAsFactors = FALSE,
                   data.table = FALSE)

write.xlsx(data.frame(st3), file="output/Supplementary_Tables_AntihypertensivesMR.xlsx", 
           sheetName="ST3",append=TRUE, row.names = FALSE, showNA = FALSE)

# Supplementary Table 4 =======================================================

mr_sbp_ad <- fread("data/mr_sbp_ad.csv",
                   stringsAsFactors = FALSE,
                   data.table = FALSE)

tmp1 <- mr_sbp_ad[!is.na(mr_sbp_ad$beta.mr),
                  c("leveli","methodtype","beta.mr","lci.mr","uci.mr","se.mr","pval.mr")]

tmp1 <- reshape(tmp1, idvar = "leveli", timevar = "methodtype", direction = "wide")

tmp2 <- mr_sbp_ad[mr_sbp_ad$methodtype=="main",c(1:3,9,13:20)]

st4 <- merge(tmp2,tmp1,all.x = TRUE)

st4 <- st4[,c(12,1,5,4,3,2,6:11,13:22)]

write.xlsx(data.frame(st4), file="output/Supplementary_Tables_AntihypertensivesMR.xlsx", 
           sheetName="ST4",append=TRUE, row.names = FALSE, showNA = FALSE)

# Supplementary Table 5 =======================================================

st5 <- unique(mr_sbp_ad[!is.na(mr_sbp_ad$egger_intercept) & mr_sbp_ad$nsnp>=10,
                        c("leveli","egger_intercept","egger_se","egger_pval","nsnp")])

colnames(st5) <- c("analysis",colnames(st4)[2:5])

write.xlsx(data.frame(st5), file="output/Supplementary_Tables_AntihypertensivesMR.xlsx", 
           sheetName="ST5",append=TRUE, row.names = FALSE, showNA = FALSE)
