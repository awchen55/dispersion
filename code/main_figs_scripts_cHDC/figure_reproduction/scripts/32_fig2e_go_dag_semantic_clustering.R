#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(GO.db)
  library(AnnotationDbi)
  library(GOSemSim)
  library(org.Hs.eg.db)
})

# ------------------------------- CONSTANTS -----------------------------------
repo_dir <- "/project/gilad/brendan/dispersion/pilot/cHDC_data"
analysis_dir <- file.path(repo_dir, "figure_reproduction")
fig2_dir <- file.path(analysis_dir, "outputs", "fig2_replication")

# rrvgo's own default reduction threshold (similarity >= 0.7 => same
# cluster); expressed here as a cutree height since we cluster on distance
# (1 - similarity), not similarity directly.
sim_threshold <- 0.7
cut_height <- 1 - sim_threshold

# ------------------------------- LOAD -----------------------------------------
csv_path <- file.path(fig2_dir, "fig2_part3_card_vs_foregut_retained.csv")
if (!file.exists(csv_path)) stop("Missing ", csv_path, " -- run 04_replicate_fig2_panels.R first")
res_df <- read.csv(csv_path)
message(sprintf("[load] %d significant terms", nrow(res_df)))

# ------------------------------- MAP MSigDB NAMES -> REAL GO IDS -----------------
bp_terms <- AnnotationDbi::select(GO.db, keys = keys(GO.db, keytype = "GOID"),
                                   columns = c("TERM", "ONTOLOGY"), keytype = "GOID")
bp_terms <- bp_terms[bp_terms$ONTOLOGY == "BP" & !is.na(bp_terms$TERM), ]
norm_go_term <- toupper(gsub("-", " ", bp_terms$TERM))
stopifnot("GO.db BP term names are not unique after hyphen normalization -- mapping would be ambiguous" =
            !any(duplicated(norm_go_term)))
name_to_goid <- setNames(bp_terms$GOID, norm_go_term)

msigdb_norm_name <- toupper(gsub("_", " ", sub("^GOBP_", "", res_df$ID)))
res_df$go_id <- unname(name_to_goid[msigdb_norm_name])
n_unmapped <- sum(is.na(res_df$go_id))
if (n_unmapped > 0) {
  stop(n_unmapped, " of ", nrow(res_df), " GOBP terms could not be mapped to a real GO ID: ",
       paste(res_df$ID[is.na(res_df$go_id)], collapse = ", "))
}
message(sprintf("[map] %d/%d terms mapped to real GO IDs (0 unmapped)", nrow(res_df), nrow(res_df)))

# ------------------------------- SEMANTIC SIMILARITY -----------------------------
# computeIC = FALSE: Wang's method is pure DAG-topology (no gene-frequency-
# based information content needed), so skip that (slower) computation.
semdata <- godata(OrgDb = org.Hs.eg.db, ont = "BP", computeIC = FALSE)
sim_mat <- mgoSim(res_df$go_id, res_df$go_id, semData = semdata, measure = "Wang", combine = NULL)
rownames(sim_mat) <- res_df$Description
colnames(sim_mat) <- res_df$Description
message(sprintf("[semsim] Wang similarity matrix: %d x %d", nrow(sim_mat), ncol(sim_mat)))

# ------------------------------- CLUSTER -----------------------------------------
dist_mat <- 1 - sim_mat
hc <- hclust(as.dist(dist_mat), method = "average")
clusters <- cutree(hc, h = cut_height)

cluster_df <- res_df %>%
  dplyr::mutate(cluster = clusters[Description]) %>%
  dplyr::arrange(cluster, p.adjust)

representatives <- cluster_df %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_min(p.adjust, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

n_clusters <- length(unique(clusters))
message(sprintf("[cluster] %d terms -> %d clusters (similarity threshold %.2f) -> %d representatives",
                nrow(res_df), n_clusters, sim_threshold, nrow(representatives)))

cluster_sizes <- cluster_df %>% dplyr::count(cluster, sort = TRUE)
message("[cluster] size distribution: ", paste(cluster_sizes$n, collapse = ", "))

# Flag which representatives are heart/cardiac-related -- purely informational.
representatives$is_heart_related <- grepl("HEART|CARDIA", toupper(representatives$Description))
message(sprintf("[cluster] %d/%d representatives are heart/cardiac-related",
                sum(representatives$is_heart_related), nrow(representatives)))

write.csv(cluster_df, file.path(fig2_dir, "fig2e_godag_semantic_clusters.csv"), row.names = FALSE)
write.csv(representatives, file.path(fig2_dir, "fig2e_godag_semantic_representatives.csv"), row.names = FALSE)

message("[done] wrote fig2e_godag_semantic_{clusters,representatives}.csv to ", fig2_dir)

# ------------------------------- REPORTABLE VALUES ------------------------------
write.csv(
  data.frame(
    key = c("fig2e_godag_n_significant_terms", "fig2e_godag_similarity_threshold",
            "fig2e_godag_n_clusters", "fig2e_godag_n_heart_representatives"),
    value = c(nrow(res_df), sim_threshold, n_clusters, sum(representatives$is_heart_related)),
    provenance = "figure_reproduction/scripts/32_fig2e_go_dag_semantic_clustering.R"
  ),
  file.path(analysis_dir, "outputs", "reportable_values_32.csv"),
  row.names = FALSE
)
