Analysis-ready datasets (written by BPD_meta_analysis_full.R, section 7b)

<GSE>_analysis_ready.rds   one per cohort (5 blood discovery, 2 non-blood, 2 external/supplementary). A list with:
    cohort, role, platform, gpl, tissue, sampling, bpd_definition
    expression  genes x samples matrix exactly as analysed
                (microarray: log2, probes collapsed to genes; RNA-seq: TMM-normalised log2 CPM)
                first-week samples only; GSE302178 preterm infants only
    pheno       sample, bpd (original definition), bpd_moderate_severe, severity, ga_weeks, birth_weight_g
                (GSE302178: chorioamnionitis instead of GA/BW)
    scores      per-sample BPD tracking index, six program scores and the 21-gene signature score (final model)

longitudinal_GSE32472_GSE220135.rds
    expression  all time points for both cohorts (genes x samples)
    annotation  infant, time point (0 = week 1/cord, 1 = D14, 2 = D28), severity, GA, BW
    scores      program and index scores used in the mixed models and trajectory classes

GSE220135_cell_fractions.rds   deposited deconvolution estimates (B, CD4 T, CD8 T, NK, granulocytes, monocytes, nRBC)

model_and_meta_analysis.rds    21-gene signature, GO anchors, final and leave-one-cohort-out core genes,
                               cohort-level effect sizes (blood, GA-adjusted, 7 cohorts), meta-analyses,
                               discovery and held-out AUCs, GO pathway results

BPD_all_analysis_data.rds      everything above in one list (written by the script; omitted here to halve the download)

Example in R:
    d <- readRDS("GSE32472_analysis_ready.rds"); dim(d$expression); head(d$pheno)
