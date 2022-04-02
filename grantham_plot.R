# Load required packages 
library(grantham)
library(tidyverse)
library(seqinr)

# Load protein alignment sequences into R
FAR5_SUB <- 'MELNCVQFLQNKTILVTGATGFLAKVFVEKILRVQPNVKKLYLLVRASDNEAATKRLHMEIFEKDLFKVLRENLGDEKLNTLLYEKVVSVPGDIATDQLGMNDSDLRERMQKEIDIVVNVAATTNFDERYDVGLGINTFGALNVLNFAKKCVKVELLLHVSTAYVCGEKPGLIPEKPFIMEDIRNENGLQLDINLERELMKQRLKELNEQDCSEEDITLSMKELGMERAKLHGWPNTYVFTKSMGEMLLGNHKENLPLVIIRPTMITSTLSEPFPGWIEGLRTVDSVIIAYGKGVLKCFLVDVNSVCDMIPVDMVANAMITAAAKHAGCSRVHMVYHVGSSHQNPVTFGEIHEIAVRYFTKNPLRSRNGSLITVSKVRFISTMALFSLYMTLRYKLPLQMLKLIDIIYPWRNGDKYGDKNRKIELVMRLVELYEPYVLFKGIFEDRNTKSLCANQKEEEIKNTEKLMFDFDPKGINWGDYLTNIHISGLVTHVLKK'
FAR5_ZEP <- 'MELNCVQFLQNKTILVTGATGFLAKVFVEKILRVQPNVKKLYLLVRASDNEAATKRLHMEIFEKDLFKVLRENLGDEKLNTLLYEKVVSVPGDIATDQLGMNDSDLRERMQKEIDIVVNVAATTNFDERYDVGLGINTFGALNVLNFAKKCVKVELLLHVSTAYVCGEKPGLIPEKPFIMEDIRNENGLQLDINLERELMKQRLKELNEQDCSEEEITLSMKELGMERAKLHGWPNTYVFTKSMGEMLLGNHKENLPLVIIRPTMITSTLSEPFPGWIEGLRTVDSVIIAYGKGVLKCFLVDVNSVCDMIPVDMVANAMITAAAKHAGGSGVHMVYHVGSSHQNPVTFGEIHEIAVRYFTKNPLRSRNGSLITVSKLRFISTMSLFSLYMTLRYKLPLQMLKLIDIIYPWRNGDKYGDKNRKIELAMRLVELYEPYVLFKGIFDDRNTKSLCANQKKEEIKNTEKLMFDFDPKGIKWGDYLTNIHISGLITHVLKK'
MAP18_SUB <- 'MGYWKSKVVPRMKRLFEKSPAKKEVVEEEKPREVEVVEEVVVKTEEPAKVEETKPEEIATGEKEIEIVEEKKEEAKPVEVPVPVAAEEKKLAVEEEKKTAPVEEKKLAVEEEKKPAVEEKKPVVEEKKKLLPRFRWLKLL*PRLPKLRWLKLRQRLRKLRWLRHKRLEFSSWYIFLKKY*----'
MAP18_ZEP <- 'MGYWKSKVVPRMKKLFEKNSAKKEVVEEEKPREVEVVEEVVVKTEEPAKVEETKPEEIATGEKEIEIVEEKKEEAKPVEVPVPV-AEEKKLAVEEEKKTAPVEEKKPAVEEEKKPAVEEKKPVVEEKKEVVAAVPVAETPSTKAPETPVVETPAKAPETPVAAPQKA----*IFFMVHFSKKIL'

# Split protein alignment sequences into a list of individual amino acids 
split_function <- function(x){
  strsplit(x, split='')
}

FAR5_SUB_split <- split_function(FAR5_SUB)
FAR5_ZEP_split <- split_function(FAR5_ZEP)
MAP18_SUB_split <- split_function(MAP18_SUB)
MAP18_ZEP_split <- split_function(MAP18_ZEP)

# Convert single letter amino acid code into triplet amino acid code (needed for grantham package)
triplet_code <- function(x){
    aaa(x)
  }

FAR5_SUB_AAA <- lapply(FAR5_SUB_split, triplet_code)
FAR5_ZEP_AAA <- lapply(FAR5_ZEP_split, triplet_code)
MAP18_SUB_AAA <- lapply(MAP18_SUB_split, triplet_code)
MAP18_ZEP_AAA <- lapply(MAP18_ZEP_split, triplet_code)

# Create a tibble for each protein 
FAR5 <- map2_df(FAR5_SUB_AAA, FAR5_ZEP_AAA, ~tibble(SUB=.x, ZEP=.y))
MAP18 <- map2_df(MAP18_SUB_AAA, MAP18_ZEP_AAA, ~tibble(SUB=.x, ZEP=.y))

# Add a position column to each protein's tibble 
FAR5$POS <- 1:496
MAP18$POS <- 1:184

# Add a substitution identification column to each tibble 
FAR5_substitutions <- paste(FAR5$SUB, FAR5$POS, FAR5$ZEP, sep='')
MAP18_substitutions <- paste(MAP18$SUB, MAP18$POS, MAP18$ZEP, sep='')

FAR5$Substitution_ID <- FAR5_substitutions
MAP18$Substitution_ID <- MAP18_substitutions

# Exclude unwanted information from MAP18 tibble for use downstream 
MAP18_new <- MAP18 %>%
  na.omit() %>%
  filter(!SUB=='Stp') %>%
  filter(!ZEP=='Stp')

# Calculate grantham distances for each protein and add this information to each tibble
FAR5_grantham <- grantham_distance(x=FAR5$SUB, y=FAR5$ZEP)
MAP18_grantham <- grantham_distance(x=MAP18_new$SUB, y=MAP18_new$ZEP)

FAR5$Grantham_score <- FAR5_grantham$d
MAP18_new$Grantham_score <- MAP18_grantham$d

# Create Grantham score plots for each protein 
FAR5_plot <- ggplot(data=FAR5, aes(x=POS, y=Grantham_score)) + geom_point(stat='identity') + ggtitle('Grantham scores of substitutions between \n SUB and ZEP for FAR5') + ylab('Grantham score')
FAR5_plot

MAP18_plot <- ggplot(data=MAP18_new, aes(x=POS, y=Grantham_score)) + geom_point(stat='identity') + ggtitle('Grantham scores of substitutions between \n SUB and ZEP for MAP18') + ylab('Grantham score')
MAP18_plot

# Filter out Grantham scores of 0 
FAR5_filter <- FAR5 %>%
  filter(Grantham_score>0)

MAP18_filter <- MAP18_new %>%
  filter(Grantham_score>0)

# Write table out to csv 
write_csv(FAR5_filter, '/Users/baileyfrancis/OneDrive/Documents/UoN/GroupProject2/FAR5_Grantham.csv')
write_csv(MAP18_filter, '/Users/baileyfrancis/OneDrive/Documents/UoN/GroupProject2/MAP18_Grantham.csv')
