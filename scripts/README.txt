BPD transcriptomic analyses
==========================

Run
---
- Run from this folder: Rscript run_all.R
- Place the eight analysis-ready RDS files below in ./data before running.
- Generated results are written to ./out; figures are written to ./figs.
- Reported analysis environment: R 4.4.3.
- Required packages: dplyr, tidyr, limma, pROC, nlme, msigdbr, data.table, ggplot2 and patchwork.
- The supplied analysis specifies MSigDB 7.5.1; confirm the installed msigdbr package supplies that release.
- Install Liberation Sans for the specified figure font and use an R build with Cairo support for PDF export.

Required data files
-------------------
- data/GSE32472_analysis_ready.rds
- data/GSE225881_analysis_ready.rds
- data/GSE125873_analysis_ready.rds
- data/GSE220135_analysis_ready.rds
- data/GSE302178_analysis_ready.rds
- data/longitudinal_GSE32472_GSE220135.rds
- data/GSE220135_cell_fractions.rds
- data/model_and_meta_analysis.rds

Script map
----------
- 00_setup.R: shared data loading, statistical functions and program scoring.
- 06_signature_rebuild.R and 06c_adjusted.R: signature reconstruction and GA-adjusted comparison.
- 07_gsva_pathways.R: pathway scoring with the supplied custom GSVA implementation.
- 09_programs_final.R and 09b_programs4.R: pathway clustering and four-program assignment.
- 10_heldout_programs.R: leave-one-cohort-out program reconstruction.
- 11_fig3.R and 13_fig4.R: program validation, GA matching, composition and longitudinal analyses.
- 14_fig5_analysis.R: top-gene prediction, gain beyond GA and published-signature audit.
- 15_timecourse.R: longitudinal prediction and exploratory response/composition adjustment.
- 16_spon1.R through 18_nt_spon1.R: SPON1, neutrophil-T-cell scores and combined models.
- 19_bigstory.R and 20_naiveT.R: GA associations, restricted-GA comparisons and naive T-cell exploration.
- 21_v7_partA.R and 22_v7_partB.R: protocol v7 Part A and Part B analyses.
- 23_newfigs_data.R and 24_numbers_external.R: revised figure data and external-cohort summaries.
- 27_axis_pca.R and 28_axis_null_sensitivity.R: updated signature PCA, random-gene nulls and sensitivity analyses.
- nfig1.R through nfig5.R: revised main Figures 1-5, saved with the NewFigure prefix.
- fig1.R, fig2.R and fig5.R: earlier figure builds retained in the supplied run order.
- 25_supp_figs.R: supplementary figures.
- 26_supp_tables_export.R: supplementary CSV exports to out/supp_tables.
- graphical_abstract.R: graphical abstract.
- figtheme.R: shared figure colours, typography and export helpers.
- source_section.R: named-section loading used by the updated PCA and sensitivity scripts.
- run_all.R: sequential execution in separate Rscript processes with RELAX=TRUE for training folds.

Implementation notes
--------------------
- Custom GSVA scores have not been checked here for equivalence with a Bioconductor GSVA release.
- Protocol documents, manuscript builders, analysis datasets and generated results are not included.
- REVIEW_NOTES.txt identifies existing analysis concerns without changing model definitions or results.
