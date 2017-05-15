library(data.table); library(reshape2); library(ggplot2); library(plotly); library(shiny)

### get LUAD samples only
cells <- read.delim("Cell_Lines_Details.txt")
cells <- cells[which(cells$Cancer.Type..matching.TCGA.label. == "LUAD"),]

### Get drugs
drugs <- read.delim("Screened_Compounds.txt")

### Get the responses of dabrafenib and cell lines of interest
response <- read.delim("response.txt")

### Read in the WES and look at only KRAS mutant samples
wes <- read.delim("WES_KRAS.txt")

### Find KRAS LUAD and merge with drug response
response$KRAS_status <- response$COSMIC_ID %in% wes$COSMIC_ID

### Merge drug names for targets
response <- merge(response, drugs, by.x = "DRUG_ID", by.y = "DRUG.ID")

### Remove extraneous dataframes
rm(cells); rm(drugs); rm(wes)

### reshape the data and calculate pvalues
df <- NULL
for(i in unique(response$DRUG_ID)){
     mut <- response[response$KRAS_status == "TRUE" & response$DRUG_ID == i,]
     wt <- response[response$KRAS_status == "FALSE" & response$DRUG_ID == i,]
     
     mut_mn <- mean(mut$AUC)
     wt_mn <- mean(wt$AUC)
     
     p_val <- wilcox.test(mut$AUC, wt$AUC)$p.val
     
     df <- rbind(df,
                 data.frame(DRUG_ID = i,
                            drug_common = unique(mut$DRUG.NAME),
                            drug_target = unique(mut$TARGET),
                            drug_pathway = unique(mut$TARGET.PATHWAY),
                            p_value = p_val,
                            neg_log_pval = -log10(p_val),
                            mut_mean = mut_mn,
                            wt_mean = wt_mn,
                            log_2_fc = log2(mut_mn / wt_mn),
                            abs_diff = abs(mut_mn - wt_mn)
                 ))
}

df <- df[order(df$drug_common),]


shinyUI(fluidPage(
     titlePanel("GDSC: LUAD and KRAS case study"),
     sidebarLayout(
          sidebarPanel(
               selectInput("drug", "Pick a Drug:",
                           choices = unique(df$drug_common))
          ),
     mainPanel(
          splitLayout(cellWidths = c("70%","30%"),
                      plotlyOutput("interactive_plot"),
                      plotOutput("boxplot")))
     )
     )
)
               
               
               
               
               
