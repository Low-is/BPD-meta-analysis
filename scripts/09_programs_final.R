# Group conserved pathways by the overlap of direction-matched signature genes.
# Cluster pathways with >= 3 leading genes by Jaccard distance and average linkage at cut 0.9.
# Retain clusters with >= 3 pathways and >= 8 genes; rank core genes by P and retain up to 30.
source("./00_setup.R")
load(file.path(OUT, "pathways_gsva.RData"))
load(file.path(OUT, "signature_rebuild.RData"))
cons <- PM |> filter(conserved) |> mutate(dir = ifelse(es < 0, -1, 1))
cons$lead <- lapply(seq_len(nrow(cons)),
  function(i) SIGR$gene[SIGR$gene %in% GOL[[cons$set[i]]] & sign(SIGR$es) == cons$dir[i]])
cons$n_lead <- lengths(cons$lead)
cons$cluster <- NA_character_
for (d in c(-1, 1)) {
  i <- which(cons$dir == d & cons$n_lead >= 3)
  L <- cons$lead[i]
  n <- length(L)
  J <- matrix(0, n, n)
  for (a in 1:n) for (b in 1:n) J[a,
    b] <- length(intersect(L[[a]], L[[b]])) / length(union(L[[a]], L[[b]]))
  cons$cluster[i] <- paste0(ifelse(d < 0, "L", "H"), cutree(hclust(as.dist(1 - J), "average"), h = .9))
}
PR <- cons |> filter(!is.na(cluster)) |> group_by(cluster) |>
  summarise(n_path = n(), dir = dir[1], genes = list(unique(unlist(lead))),
    top = paste(head(sub("GOBP_", "", set[order(p)]), 3), collapse = " | ")) |>
  mutate(n_lead = lengths(genes)) |> filter(n_path >= 3, n_lead >= 8) |> arrange(dir, desc(n_lead))
PR$core <- lapply(PR$genes, function(g) SIGR |> filter(gene %in% g) |> arrange(p) |> head(30))
print(PR |> select(cluster, dir, n_path, n_lead, top), width = 250)
for (i in seq_len(nrow(PR))) cat(PR$cluster[i], ":", paste(head(PR$core[[i]]$gene, 10), collapse = ", "),
  "\n")
# overlap of core genes between programs
ov <- outer(seq_len(nrow(PR)), seq_len(nrow(PR)),
  Vectorize(function(a, b) length(intersect(PR$core[[a]]$gene, PR$core[[b]]$gene))))
dimnames(ov) <- list(PR$cluster, PR$cluster)
print(ov)
save(cons, PR, file = file.path(OUT, "programs_final.RData"))
