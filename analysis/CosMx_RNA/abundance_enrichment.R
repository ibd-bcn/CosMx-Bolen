library(ggplot2)
library(Seurat)
library(readr)
library(dplyr)

todas <-   readRDS("analysis/CosMx_RNA/Objects/seurats_all_norm.RDS")

## IBD vs NHC ------------------------------------------------------------------
df_final <- data.frame(matrix(ncol = 4, nrow = 0))
cnames <- c("Cluster", "Response", "stdres", "P.value")
colnames(df_final) <- cnames

long.data <-
  todas@meta.data[todas@meta.data$tissue %in% c("IBD", "NHC"), ]


freq.table <- long.data %>%
  dplyr::group_by(tissue, patient, refined) %>%
  dplyr::tally()

chis.table <- table(long.data$tissue, long.data$refined)

# loop for each row:
for (i in 1:ncol(chis.table)) {
  Chi2 <- chisq.test(chis.table[, i])

  wilx <- wilcox.test(n ~ tissue, data = freq.table[freq.table$refined == colnames(chis.table)[i],])


  stdres <- Chi2$stdres
  pval <- Chi2$p.value
  fc <- chis.table[1, i] / chis.table[2, i]

  df <- data.frame(
    row.names = colnames(chis.table)[i],
    "Cluster" = colnames(chis.table)[i],
    "e.score" = Chi2$stdres[grep('IBD', names(Chi2$stdres))],
    "Comparison" = "IBD",
    "Wlx_p.value" = wilx$p.value,
    "Chisq_p.value" = Chi2$p.value,
    'FC' = fc
  )

  df_final <- rbind(df_final, df)

}



df_final$fdr_wlx <- p.adjust(df_final$Wlx_p.value, method = 'BY')

df_final$ast_wlx <- ifelse(
  df_final$FC > 1.5 | df_final$FC < 0.66,
  yes = ifelse(
    df_final$fdr_wlx < 0.05,
    no = '',
    yes = ifelse(
      df_final$fdr_wlx < 0.01,
      no = '*',
      yes = ifelse(df_final$fdr_wlx < 0.001,
                   no = '**',
                   yes = '***')
    )
  ),
  no = ''
)

df_final$fdr_Chisq <-
  p.adjust(df_final$Chisq_p.value, method = 'BY')
df_final$ast_Chisq <- ifelse(
  df_final$FC > 1.5 | df_final$FC < 0.66,
  yes = ifelse(
    df_final$fdr_Chisq < 0.05,
    no = '',
    yes = ifelse(
      df_final$fdr_Chisq < 0.01,
      no = '*',
      yes = ifelse(df_final$fdr_Chisq < 0.001,
                   no = '**',
                   yes = '***')
    )
  ),
  no = ''
)

df_final$log10pval_wlx <- -log10(df_final$fdr_wlx)
df_final$log10pval_Chisq <- -log10(df_final$fdr_Chisq)

## PD vs NHC -------------------------------------------------------------------

df_final1 <- data.frame(matrix(ncol = 4, nrow = 0))
cnames <- c("Cluster", "Response", "stdres", "P.value")
colnames(df_final1) <- cnames

long.data <-
  todas@meta.data[todas@meta.data$tissue %in% c("PD", "NHC"), ]


freq.table <- long.data %>%
  dplyr::group_by(tissue, patient, refined) %>%
  dplyr::tally()

chis.table <- table(long.data$tissue, long.data$refined)

# loop for each row:
for (i in 1:ncol(chis.table)) {
  Chi2 <- chisq.test(chis.table[, i])

  wilx <- wilcox.test(n ~ tissue, data = freq.table[freq.table$refined == colnames(chis.table)[i],])


  stdres <- Chi2$stdres
  pval <- Chi2$p.value
  fc <- chis.table[1, i] / chis.table[2, i]

  df <- data.frame(
    row.names = colnames(chis.table)[i],
    "Cluster" = colnames(chis.table)[i],
    "e.score" = Chi2$stdres[grep('PD', names(Chi2$stdres))],
    "Comparison" = "PD",
    "Wlx_p.value" = wilx$p.value,
    "Chisq_p.value" = Chi2$p.value,
    'FC' = fc
  )

  df_final1 <- rbind(df_final1, df)

}



df_final1$fdr_wlx <- p.adjust(df_final1$Wlx_p.value, method = 'BY')

df_final1$ast_wlx <- ifelse(
  df_final1$FC > 1.5 | df_final1$FC < 0.66,
  yes = ifelse(
    df_final1$fdr_wlx < 0.05,
    no = '',
    yes = ifelse(
      df_final1$fdr_wlx < 0.01,
      no = '*',
      yes = ifelse(df_final1$fdr_wlx < 0.001,
                   no = '**',
                   yes = '***')
    )
  ),
  no = ''
)

df_final1$fdr_Chisq <-
  p.adjust(df_final1$Chisq_p.value, method = 'BY')
df_final1$ast_Chisq <-
  ifelse(
    df_final1$FC > 1.5 | df_final1$FC < 0.66,
    yes = ifelse(
      df_final1$fdr_Chisq < 0.05,
      no = '',
      yes = ifelse(
        df_final1$fdr_Chisq < 0.01,
        no = '*',
        yes = ifelse(df_final1$fdr_Chisq < 0.001,
                     no = '**',
                     yes = '***')
      )
    ),
    no = ''
  )

df_final1$log10pval_wlx <- -log10(df_final1$fdr_wlx)
df_final1$log10pval_Chisq <- -log10(df_final1$fdr_Chisq)


## Prepare for Plots -----------------------------------------------------------
df_final <- rbind(df_final, df_final1)
df_final$Cluster <- factor(df_final$Cluster)
df_final$Comparison <- factor(df_final$Comparison)

write.csv("analysis/CosMx_RNA/Results/enrichment.csv")

