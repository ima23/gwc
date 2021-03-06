studies = list()
studies$gnomad = list("study_name" = "gnomad", "vars" = "~/scratch/gnomAD/gnomad.genomes.r2.0.1.sites.%s.txt", 'null_model' = "~/software/selection2/models/obs_exp_lm.gnomAD_NFE.30x_cov.RData", "pop_size" = 7509, "AC_column_name" = 'AC_NFE', "gnomad_cov_correction" = "~/software/selection2/models/gnomad_cov_correction_loess.RData", "pass_flags" = "PASS", "neutral_lm" = "~/software/selection2/models/obs_exp_lm.gnomAD.phylop.RData", maps_lm = "~/software/selection2/models/maps_lm.gnomAD.phylop.RData")
studies$BRIDGE = list("study_name" = "BRIDGE", "vars" = "~/scratch/BRIDGE/chr%s.tsv", 'null_model' = "~/software/selection2/models/obs_exp_lm.BRIDGE.30x_cov.RData", "pop_size" = 13049, "AC_column_name" = 'AC', "bridge_cov_correction" = "~/software/selection2/models/bridge_cov_correction_loess.RData", "pass_flags" = "PASS", neutral_lm = "~/software/selection2/models/obs_exp_lm.BRIDGE.phylop.RData", maps_lm = "~/software/selection2/models/maps_lm.BRIDGE.phylop.RData", neutral_lm_no_methyl = "~/software/selection2/models/obs_exp_lm.no_methyl_correction.BRIDGE.phylop.RData")

mut_model_covariates = list()
mut_model_covariates$low_complexity_regions = list("mark_name" = "low_complexity_regions", "file" = "~/scratch/MUTMODEL/LCR-hs37d5.bed.gz", method = "overlap")
mut_model_covariates$replication_timing_LCLs = list("mark_name" = "replication_timing_Koren_LCLs", "file" = "~/scratch/MUTMODEL/Koren_LCL_rep_timing.txt", method = "nearest")
mut_model_covariates$replication_timing_ESC = list("mark_name" = "replication_timing_DingKoren_ESCs", "file" = "~/scratch/MUTMODEL/peaks/ESC.rep_timing.koren_ding.averaged.sorted.bed.gz", 'method' = 'nearest')
mut_model_covariates$recombination_rate_kong_female = list("mark_name" = "recombination_rate_kong_female", "file" = "~/scratch/MUTMODEL/kong_2010_female_rmap.bed", method = "nearest")
mut_model_covariates$recombination_rate_kong_male = list("mark_name" = "recombination_rate_kong_male", "file" = "~/scratch/MUTMODEL/kong_2010_male_rmap.bed", method = "nearest")
mut_model_covariates$recombination_rate_1kg = list("mark_name" = "recombination_rate_1000G_phase3", file = "~/scratch/MUTMODEL/peaks/IMPUTE2_1000G_Phase3_recombination_map.sorted.bed.gz", method = 'nearest')
mut_model_covariates$CTCF_binding_site = list("mark_name" = "CTCF_binding_site", file = "~/scratch/MUTMODEL/homo_sapiens.GRCh37_liftover.Regulatory_Build.CTCF_binding_sites.20161111.sorted.bed", method = "overlap")
mut_model_covariates$ovary_DNase = list("mark_name" = "ovary_DNase", file =	"~/scratch/MUTMODEL/peaks/E097-DNase.hotspot.fdr0.01.broad.sorted.bed.gz", method = "overlap")
mut_model_covariates$hSSC_ATAC = list("mark_name" = "hSSC_ATAC", file = "~/scratch/MUTMODEL/peaks/hSSC_atac.sorted.bed.gz", method = 'overlap')
mut_model_covariates$hESC_ATAC = list("mark_name" = "hESC_ATAC", file =	"~/scratch/MUTMODEL/peaks/hESC_atac.sorted.bed.gz", method = 'overlap')
mut_model_covariates$ovary_H3K27ac = list("mark_name" = "ovary_H3K27ac", file = "~/scratch/MUTMODEL/peaks/E097-H3K27ac.sorted.broadPeak.gz", method = "overlap")
mut_model_covariates$ovary_H3K27me3 = list("mark_name" = "ovary_H3K27me3", file =	"~/scratch/MUTMODEL/peaks/E097-H3K27me3.sorted.broadPeak.gz", method = "overlap")
mut_model_covariates$ovary_H3K9me3 = list("mark_name" = "ovary_H3K9me3", file =	"~/scratch/MUTMODEL/peaks/E097-H3K9me3.sorted.broadPeak.gz", method = "overlap")
mut_model_covariates$ovary_H3K4me3 = list("mark_name" = "ovary_H3K4me3", file =	"~/scratch/MUTMODEL/peaks/E097-H3K4me3.sorted.broadPeak.gz", method = "overlap")
mut_model_covariates$ovary_H3K4me1 = list("mark_name" = "ovary_H3K4me1", file = "~/scratch/MUTMODEL/peaks/E097-H3K4me1.sorted.broadPeak.gz", method = "overlap")
mut_model_covariates$ovary_H3K36me3 = list("mark_name" = "ovary_H3K36me3", file = "~/scratch/MUTMODEL/peaks/E097-H3K36me3.sorted.broadPeak.gz", method = "overlap")
mut_model_covariates$hESC_H3K27me3 = list("mark_name" = "hESC_H3K27me3", file =	"~/scratch/MUTMODEL/peaks/E001-H3K27me3.sorted.broadPeak.gz", method = "overlap")
mut_model_covariates$hESC_H3K9me3 = list("mark_name" = "hESC_H3K9me3", file =	"~/scratch/MUTMODEL/peaks/E001-H3K9me3.sorted.broadPeak.gz", method = "overlap")
mut_model_covariates$hESC_H3K4me3 = list("mark_name" = "hESC_H3K4me3", file =	"~/scratch/MUTMODEL/peaks/E001-H3K4me3.sorted.broadPeak.gz", method = "overlap")
mut_model_covariates$hESC_H3K4me1 = list("mark_name" = "hESC_H3K4me1", file = "~/scratch/MUTMODEL/peaks/E001-H3K4me1.sorted.broadPeak.gz", method = "overlap")
mut_model_covariates$hESC_H3K36me3 = list("mark_name" = "hESC_H3K36me3", file = "~/scratch/MUTMODEL/peaks/E001-H3K36me3.sorted.broadPeak.gz", method = "overlap")
mut_model_covariates$hESC_H3K9ac = list("mark_name" = "hESC_H3K9ac", file = "~/scratch/MUTMODEL/peaks/E001-H3K9ac.sorted.broadPeak.gz", method = "overlap")
mut_model_covariates$ovary_GTex = list("mark_name" = "ovary_GTex", file = "~/scratch/MUTMODEL/gencode_v19_transcripts_ovary_expression_GTex.bed", method = "overlap_value")
mut_model_covariates$testis_GTex = list("mark_name" = "testis_GTex", file = "~/scratch/MUTMODEL/gencode_v19_transcripts_testis_expression_GTex.bed", method = "overlap_value")
mut_model_covariates$gencode_v19_protein_coding = list("mark_name" = "gencode_v19_protein_coding", file = "~/scratch/MUTMODEL/gencode.v19.CDS.protein_coding.bed", method = "overlap")
