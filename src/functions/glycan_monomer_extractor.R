glycan_monomer_extractor <- function (x) {
  tmp_1 <- strsplit(gsub("\\(N-Glycosylation\\)|\\(Oxidation\\)", "*", x), "\\{")[[1]][2]
  tmp_2 <- gsub("[[:punct:]]","",tmp_1)
  tmp_3 <- as.numeric(c(
    str_match(tmp_2, "Hex(.)")[2],
    str_match(tmp_2, "HexNAc(.)")[2],
    str_match(tmp_2, "Neu5Ac(.)")[2],
    str_match(tmp_2, "Fuc(.)")[2]))
  names(tmp_3) <- c("Hex","HexNAc","Neu5Ac","Fuc")
  return(tmp_3)
}