# Export supplementary tables (CSV) for assembly into the supplementary workbook
source("./00_setup.R")
for (f in c("signature_rebuild", "programs4", "fig3", "fig4", "fig5", "v7_partA", "v7_partB", "newfigs",
  "timecourse", "external_v7")) load(file.path(OUT, paste0(f, ".RData")))
TD <- "./out/supp_tables"
dir.create(TD, showWarnings = FALSE)

w <- function(d, n) write.csv(d |> mutate(across(where(is.numeric), ~ signif(.x, 4))),
  file.path(TD, paste0(n, ".csv")), row.names = FALSE, na = "")
w(SIGR |> transmute(gene, direction = ifelse(es > 0, "Higher in BPD", "Lower in BPD"), pooled_g = es, se,
  p, fdr, I2, cohorts = k, cohorts_up = n_pos, cohorts_down = n_neg) |> arrange(p), "S3_signature_1151")
cons2 <- cons |> mutate(across(where(is.list), ~ sapply(.x, function(v) paste(v, collapse = ";"))))
w(cons2, "S6_conserved_pathways")
w(bind_rows(lapply(seq_len(nrow(PR4)),
  function(i) PR4$core[[i]] |> mutate(program = PR4$program[i], .before = 1))), "S7_program_core_genes")
w(AUCT |> rename(analysis = type), "S8a_program_heldout_auc")
w(DEVT, "S8b_program_GA_matched_and_adjusted")
w(EFF, "S9a_programs_by_timepoint")
w(MM, "S9b_programs_mixed_models")
w(FI, "S9c_trajectory_classes")
SRC <- c("Luo 2024, PLoS One" = "GSE32472 (assumed; not verified)",
  "Luo 2024, Medicine" = "GSE32472 (per review)", "Fang 2024" = "GSE32472 (verified)",
  "Zou 2024" = "GSE32472 (assumed; not verified)",
  "Zhang 2025" = "GSE32472 day 28 (verified; validated on days 5 and 14)",
  "Guo 2025" = "GSE220135 (verified); validated by its authors in GSE32472",
  "Tu 2025, Front Pediatr" = "GSE32472 days 5/14 (verified); validated in GSE188944 = GSE220135 infants",
  "Sun 2026" = "GSE32472 (verified); validated in GSE108754",
  "Jia 2023" = "GSE32472 (assumed; not verified)", "Li 2024 (APP)" = "GSE32472 (assumed; not verified)",
  "Abilbayeva 2026" = "GSE32472 (verified)",
  "Ann Med 2024 (autophagy)" = "GSE32472 training set (stated in the paper)",
  "Windhorst 2023" = "GSE225881 (authors' cohort); directions from published PAM class scores",
  "Authors 2023, 5 genes" = "GSE32472 (authors' study)",
  "Authors 2025, Th2 (asthma)" = "A priori (asthma); directions from GSE32472",
  "De novo top 5 (GSE32472)" = "GSE32472 (this study)",
  "De novo top 50 (GSE32472)" = "GSE32472 (this study)")
A2 <- AUD |> mutate(derivation = SRC[signature],
  .after = signature) |> rename(genes_measured_in_test_cohorts = cov, auc_GSE32472 = auc_32,
    auc_GSE225881 = auc_81, auc_GSE125873 = auc_73, auc_GSE220135 = auc_35,
    pooled_auc = auc, auc_lo = lo, auc_hi = hi, GA_adj_effect = ga_g, GA_adj_lo = ga_lo,
    GA_adj_hi = ga_hi, GA_adj_p = ga_p, r_neutrophil_score = rN, r_Tcell_score = rT)
w(A2, "S10a_published_signature_audit")
w(SIGMAP, "S10b_signature_axes_and_balance")
w(BEYP |> rename(score = method, dAUC = d), "S11a_prediction_beyond_GA_pooled")
w(BEY, "S11b_prediction_beyond_GA_by_cohort")
w(B2, "S11c_balance_components")
w(H5, "S11d_prespecified_scores_prediction")
w(data.frame(direction = names(XPRED), t(sapply(XPRED, function(v) v[c(1, 2, 3, 4, 5)])),
  check.names = FALSE) |> setNames(c("tested_in", "dAUC", "lo", "hi", "auc_GA", "auc_GA_module")),
  "S11e_cross_cohort_module_prediction")
w(TOPN, "S11f_signature_size")
w(GZ, "S11g_gray_zone")
w(SCAT |> transmute(gene, t_GSE32472 = t32, p_GSE32472 = p32, t_GSE220135 = t35, p_GSE220135 = p35,
  naive_T_gene = naive) |> arrange(pmax(p_GSE32472, p_GSE220135)), "S12_GA_adjusted_genes_two_cohorts")
w(bind_rows(lapply(names(out),
  function(k) out[[k]] |> transmute(discovered_in = k, gene, t, p = P.Value, fdr = adj.P.Val))),
  "S13a_cross_cohort_modules")
w(XC, "S13b_cross_cohort_replication")
w(WK1, "S14a_protocol_v7_week1")
w(LONG, "S14b_protocol_v7_longitudinal")
w(B3, "S15a_declared_gene_pairs")
w(B5, "S15b_balance_vs_published_PC1")
w(B5b, "S15c_signatures_after_GA_and_balance")
w(METH |> select(id, infant, t, y, ga, Gran, CD4T, CD8T, B, NK, Mono, nRBC, mNLR, BAL),
  "S15d_methylation_vs_balance_GSE220135")
w(PRED, "S16a_timecourse_prediction")
w(GCT, "S16b_glucocorticoid_proxy")
w(ADJ, "S16c_timecourse_adjustment")
w(data.frame(score = names(EXT$auc), auc_GSE302178 = EXT$auc), "S17_external_GSE302178")
cat(length(list.files(TD)), "tables\n")
