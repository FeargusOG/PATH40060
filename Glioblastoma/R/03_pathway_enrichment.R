library(ggplot2)
ensure_package("clusterProfiler")
library(clusterProfiler)
ensure_package("org.Hs.eg.db")
library(org.Hs.eg.db)
ensure_package("enrichplot")
library(enrichplot)
ensure_package("ReactomePA")
library(ReactomePA)
ensure_package("pathview")
library(pathview)

plot_pathway_enrichment_results <- function(results, title) {
  if (!is.null(results) && nrow(results) > 0) {
    print(dotplot(results) + ggtitle(title))
    
    # Tree plot needs at least 2 results
    pw_res <- pairwise_termsim(results)
    if (nrow(pw_res) >= 2) { # Ensure at least 2 terms for clustering
      print(treeplot(pw_res) + ggtitle(title))
    } else {
      message("Treeplot cannot be generated for ", title)
    }
  } else {
    message("No significant results for ", title)
  }
}

run_pathway_enrichment <- function(deg_gene_symbols, gene_list_descr, p_thresh = 0.05, q_thresh = 0.2) {
  
  print("Gene Ontology...")
  
  ## Gene Ontology (GO) enrichment analysis
  go_results <- enrichGO(
    gene          = deg_gene_symbols,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP", # Biological Process
    pvalueCutoff = p_thresh,
    qvalueCutoff  = q_thresh
  )
  
  plot_pathway_enrichment_results(go_results, 
                                  paste("Gene Ontology Enrichment", gene_list_descr))
  
  
  print("KEGG...")
  
  ## KEGG Pathway Enrichment Analysis
  ### Map gene symbol to gene entrez ID. 
  ### This is necessary for Reactome and Keggs.
  deg_gene_entrez <- bitr(
    deg_gene_symbols,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  )
  
  ### defaults to 'hsa' for organism i.e. home sapien
  kegg_results =  enrichKEGG(
    gene          = deg_gene_entrez[,2],
    pvalueCutoff = p_thresh,
    qvalueCutoff  = q_thresh
  )
  
  plot_pathway_enrichment_results(kegg_results, 
                                  paste("Kegg Pathway Enrichment", gene_list_descr))
  
  print("Reactome...")
  
  ## Reactome Pathway Enrichment Analysis
  ### defaults to 'human' for organism
  reactome_results =  enrichPathway(
    gene          = deg_gene_entrez[,2],
    pvalueCutoff = p_thresh,
    qvalueCutoff  = q_thresh
  )
  
  plot_pathway_enrichment_results(reactome_results, 
                                  paste("Reactome Pathway Enrichment", gene_list_descr))
}
