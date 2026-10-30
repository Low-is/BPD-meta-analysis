# Build programs from the rule-based pathway clusters.
# Assign remaining pathways with >= 3 leading genes to the closest same-direction program at Jaccard >= 0.2.
# Rank the union of leading genes by P and retain up to 30 core genes.
source("./00_setup.R")
load(file.path(OUT, "pathways_gsva.RData"))
load(file.path(OUT, "signature_rebuild.RData"))
load(file.path(OUT, "programs_final.RData"))
MAP <- c(L3 = "RNA processing and translation", L9 = "RNA processing and translation",
  L7 = "RNA processing and translation",
  L4 = "RNA processing and translation", L5 = "T cell activation", H1 = "Cell cycle", H4 = "Cell cycle",
  H2 = "Neutrophil antimicrobial")
P4 <- c("RNA processing and translation", "T cell activation", "Cell cycle", "Neutrophil antimicrobial")
cons$program <- unname(MAP[cons$cluster])

pg_genes <- function(pg) unique(unlist(cons$lead[which(cons$program == pg)]))
seed <- lapply(setNames(P4, P4), pg_genes)
pdir <- sapply(P4, function(pg) cons$dir[which(cons$program == pg)][1])
cons$assigned <- FALSE
for (i in which(is.na(cons$program) & cons$n_lead >= 3)) {
  cand <- P4[pdir == cons$dir[i]]
  j <- sapply(cand,
    function(pg) length(intersect(cons$lead[[i]], seed[[pg]])) / length(union(cons$lead[[i]],
      seed[[pg]])))
  if (max(j) >= .2) {
    cons$program[i] <- cand[which.max(j)]
    cons$assigned[i] <- TRUE
  }
}
# each gene belongs to one program: the one where it is a leading gene in the most member pathways
cnt <- bind_rows(lapply(which(!is.na(cons$program)),
  function(i) tibble(gene = cons$lead[[i]], program = cons$program[i]))) |> count(gene, program) |>
  group_by(gene) |> slice_max(n, n = 1, with_ties = FALSE) |> ungroup()

pg_genes <- function(pg) cnt$gene[cnt$program == pg]
PR4 <- tibble(program = P4) |> rowwise() |> mutate(n_path = sum(cons$program == program, na.rm = TRUE),
  n_assigned = sum(cons$program == program & cons$assigned, na.rm = TRUE),
  genes = list(pg_genes(program)), n_lead = length(genes),
  core = list(SIGR |> filter(gene %in% genes) |> arrange(p) |> head(30))) |> ungroup()
print(PR4 |> select(program, n_path, n_assigned, n_lead))
for (pg in P4) {
  cat("\n", pg, ":\n  pathways:",
    paste(sub("GOBP_", "",
      cons$set[which(cons$program == pg)][order(cons$p[which(cons$program == pg)])]), collapse = "; "),
    "\n  core top 12:", paste(head(PR4$core[[which(PR4$program == pg)]]$gene, 12), collapse = ", "),
    "\n")
}
cat("\nConserved pathways not in any program:", sum(is.na(cons$program)), "\n")
save(cons, PR4, P4, file = file.path(OUT, "programs4.RData"))
