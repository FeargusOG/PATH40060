library(readr)
gbm_tcga_pub2013_cox_results <- read_csv("gbm_tcga_pub2013_cox_results.csv")

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

# Extract top genes with q-value below a threshold (e.g. 0.1)
selected_genes <- gbm_tcga_pub2013_cox_results$qvalue < 0.5
entrez_ids <- gbm_tcga_pub2013_cox_results$Entrez_Gene_Id[selected_genes]

# Run GO enrichment
ego <- enrichGO(gene          = entrez_ids,
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

ekegg <- enrichKEGG(gene         = entrez_ids,
                    organism     = 'hsa',
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

ereactome <- enrichPathway(gene         = entrez_ids,
                           organism     = "human",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05,
                           readable     = TRUE)


# View results
dotplot(ego, showCategory=10, title="GO BP Enrichment")
dotplot(ekegg, showCategory=10, title="KEGG Pathway Enrichment")
dotplot(ereactome, showCategory=10, title="Reactome Pathway Enrichment")