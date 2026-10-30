# Reproduce all results, tables and figures. Run from this folder:  Rscript run_all.R
# Expects the analysis-ready .rds files in ./data and creates ./out and ./figs.
dir.create("out", showWarnings = FALSE)
dir.create("figs", showWarnings = FALSE)
steps <- c("06_signature_rebuild.R", "06c_adjusted.R", "07_gsva_pathways.R", "09_programs_final.R",
  "09b_programs4.R",
  "10_heldout_programs.R", "11_fig3.R", "13_fig4.R", "fig1.R", "fig2.R", "14_fig5_analysis.R", "fig5.R",
  "15_timecourse.R", "16_spon1.R", "17_ntscore.R", "18_nt_spon1.R", "19_bigstory.R", "20_naiveT.R",
  "21_v7_partA.R", "22_v7_partB.R", "23_newfigs_data.R", "24_numbers_external.R", "28_axis_pca.R",
  "29_axis_null_sensitivity.R",
  "nfig1.R", "nfig2.R", "nfig3.R", "nfig4.R", "nfig5.R", "25_supp_figs.R", "26_supp_tables_export.R",
  "27_graphical_abstract.R")
Sys.setenv(RELAX = "TRUE")  # 10_heldout_programs.R: nominal criteria in training folds (see STAR Methods)
for (s in steps) {
  cat("\n====", s, "====\n")
  st <- system2("Rscript", s)
  if (st != 0) stop("failed: ", s)
}
