# pass a tsv
# and the set of elements to calculate constraint in
# and variants from unaffected parents

# uses the dddMAPS library

library(stringr)
library(plyr)
library(optparse)
source("./annotation_tools.R")
source('./mutation_null_model.R')
library(BSgenome.Hsapiens.UCSC.hg19)

option_list <- list(
  make_option("--elements", help = "Elements to annotate with constraint scores - expects 1-based coordinates."),
  make_option("--chromosome", help = "Chromosome to use (for running in parallel)"),
  make_option("--phylop_neutral_max", default = 0, help = "Cutoff for site to be considered neutrally evolving for building the model"),
  make_option("--config_file", help = "File that describes the list of studies to use. Name of the study will be used to name columns."),
  make_option("--output", help = "where to save the elements annotated with obs/exp and z scores"),
  make_option("--coverage_correction", default=F, action='store_true', help = "use coverage correction models (in config file) to make coverage correction to observed variants"),
  make_option("--bed_file", default=F, action='store_true', help = "use this option if elements is bed file. assumes first three columns are chr, start, stop"),
  make_option("--rds_file", default=F, action='store_true', help = "use this option if elements is an RDS file - used in BRIDGE regulome project")
)

args <- parse_args(OptionParser(option_list=option_list))

source(args$config_file)

# load in regions file with required columns: chr, start, stop
if (args$bed_file) {
  elements <- read.table(args$elements, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(elements)[1:3] = c("chr", "start", "end")
} else if (args$rds_file) {
  elements <- readRDS(args$elements)
} else {
  if (summary( file(args$elements) )$class == "gzfile") {
    elements <- read.table(gzfile(args$elements), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  } else {
    elements <- read.delim(args$elements)
  }
}

# make all of the in the form chr1, chr2, etc.
if (any(!grepl("^chr", elements$chr))) {
  elements$chr = paste0("chr", elements$chr)
}

if (!(grepl("chr", args$chromosome))) {
  elements = subset(elements, chr %in% paste0("chr", args$chromosome))
} else {
  elements = subset(elements, chr %in% args$chromosome)
}


if (nrow(elements) == 0){
  stop("Did not find any elements after filtering by chromosome.")
}

# add region id
if (!("region_id" %in% colnames(elements))) {
  elements$region_id = paste0(elements$chr, ":", elements$start, "-", elements$end)
}

print("Elements passed:")
print(head(elements))

# add phastcons100 score
#library(phastCons100way.UCSC.hg19)
#element_intervals = GRanges(seqnames=elements$chr, IRanges(start = elements$start, width = elements$end - elements$start + 1))
#elements$phastcons100 = scores(phastCons100way.UCSC.hg19, element_intervals)

print(sprintf("Finished annotating with sequence and phastcons score at %s", Sys.time()))

elements$filter = "PASS"


print(sprintf("Annotating elements with variant counts from each of the studies in the config file, started at %s", Sys.time()))

for (study in studies) {
  print(study)
  starting_cols = colnames(elements)  

  # load in the neutral model
  load(study$neutral_lm)  # loads as synonymous_lm

  # calculate the number of expected variants in neutral sites
  #elements$expected_neutral = predict(neutral_lm, newdata = data.frame(p_snp_phylop_lt_0 = elements$p_snp_phylop_lt_0))  

  # calculate the number of expected variants in non-neutral sites
  #elements$expected_non_neutral = predict(neutral_lm, newdata = data.frame(p_snp_phylop_lt_0 = elements$p_snp_phylop_gt_0))

  # calculate the total expected
  #elements$expected_total = predict(neutral_lm, newdata = data.frame(p_snp_phylop_lt_0 = elements$p_snp_total))

  # load in the variants to get actual observed
  vars = read.delim(sprintf(study$vars, args$chromosome), stringsAsFactors = F)
  vars$phylop = as.numeric(vars$phylop)

  if (any(!grepl("^chr", vars$chr))) {
    vars$chr = paste0("chr", vars$chr)
  }
    
  print(head(vars))

  # add observed per element
  v = filter_with_bed(vars, elements)
  print('vars pre allele count filter')
  print(head(vars))

  v = get_region_id_multi_overlap(v, elements)

  ps = ddply(v, "region_id", function(df) data.frame(number_singleton = sum(df[,study$AC_column_name] == 1), number_total = nrow(df)))

  mask = (v[,study$AC_column_name] < study$pop_size*2*0.001) & (v[,study$AC_column_name] > 0)
  v = subset(v, mask)
  print('vars post allele count filter')
  print(head(v))

  #v$neutral = v$phylop <= args$phylop_neutral_max

  print('removing all sites with phylop < -1')
  v = subset(v, phylop > -1)
  v$neutral = v$phylop > -1 & v$phylop < 1

  print('Neutral v non-neutral counts')
  print(table(v$neutral))

  v_indel = subset(v, (nchar(as.character(v$alt)) != 1) | (nchar(as.character(v$ref)) != 1))
  v = subset(v, (nchar(as.character(v$alt)) == 1) & (nchar(as.character(v$ref)) == 1)) 

  o = ddply(v, "region_id", function(df) data.frame(observed_neutral = sum(df$filter[df$neutral] %in% study$pass_flags), observed_non_neutral = sum(df$filter[!df$neutral] %in% study$pass_flags), observed_neutral_low_qual = sum(!(df$filter[df$neutral] %in% study$pass_flags)), observed_non_neutral_low_qual = sum(!(df$filter[!df$neutral] %in% study$pass_flags)), observed_CG_neutral = sum((df$filter[df$neutral] %in% study$pass_flags) & ((df$ref[df$neutral] == 'G' & df$alt[df$neutral] == 'C') | (df$ref[df$neutral] == 'C' & df$alt[df$neutral] == 'G'))), observed_CT_neutral = sum((df$filter[df$neutral] %in% study$pass_flags) & ((df$ref[df$neutral] == 'C' & df$alt[df$neutral] == 'T') | (df$ref[df$neutral] == 'G' & df$alt[df$neutral] == 'A'))), observed_CA_neutral = sum((df$filter[df$neutral] %in% study$pass_flags) & ((df$ref[df$neutral] == 'C' & df$alt[df$neutral] == 'A') | (df$ref[df$neutral] == 'G' & df$alt[df$neutral] == 'T'))), observed_TA_neutral = sum((df$filter[df$neutral] %in% study$pass_flags) & ((df$ref[df$neutral] == 'T' & df$alt[df$neutral] == 'A') | (df$ref[df$neutral] == 'A' & df$alt[df$neutral] == 'T'))),observed_TG_neutral = sum((df$filter[df$neutral] %in% study$pass_flags) & ((df$ref[df$neutral] == 'T' & df$alt[df$neutral] == 'G') | (df$ref[df$neutral] == 'A' & df$alt[df$neutral] == 'C'))), observed_TC_neutral = sum((df$filter[df$neutral] %in% study$pass_flags) & ((df$ref[df$neutral] == 'T' & df$alt[df$neutral] == 'C') | (df$ref[df$neutral] == 'A' & df$alt[df$neutral] == 'G')))))
  o_indel = ddply(v_indel, "region_id", function(df) data.frame(number_singleton = sum(df[,study$AC_column_name] == 1), number_total = nrow(df), observed_neutral_indels = sum(df$filter[df$neutral] %in% study$pass_flags), observed_non_neutral_indels = sum(df$filter[!df$neutral] %in% study$pass_flags)))

  # add the quantities relevant to singleton proportion
  elements$number_singleton = 0
  elements$number_singleton[match(ps$region_id, elements$region_id)] = ps$number_singleton

  elements$number_total = 0
  elements$number_total[match(ps$region_id, elements$region_id)] = ps$number_total

  # get number of observed variants
  elements$observed_neutral = 0
  elements$observed_neutral[match(o$region_id, elements$region_id)] = o$observed_neutral

  # add the number of observed indels
  elements$observed_neutral_indels = 0
  elements$observed_neutral_indels[match(o_indel$region_id, elements$region_id)] = o_indel$observed_neutral_indels

  elements$observed_non_neutral_indels = 0
  elements$observed_non_neutral_indels[match(o_indel$region_id, elements$region_id)] = o_indel$observed_non_neutral_indels  

  elements$observed_CA_neutral = 0
  elements$observed_CA_neutral[match(o$region_id, elements$region_id)] = o$observed_CA_neutral

  elements$observed_CG_neutral = 0  
  elements$observed_CG_neutral[match(o$region_id, elements$region_id)] = o$observed_CG_neutral

  elements$observed_CT_neutral = 0
  elements$observed_CT_neutral[match(o$region_id, elements$region_id)] = o$observed_CT_neutral

  elements$observed_TA_neutral = 0
  elements$observed_TA_neutral[match(o$region_id, elements$region_id)] = o$observed_TA_neutral

  elements$observed_TC_neutral = 0
  elements$observed_TC_neutral[match(o$region_id, elements$region_id)] = o$observed_TC_neutral

  elements$observed_TG_neutral = 0
  elements$observed_TG_neutral[match(o$region_id, elements$region_id)] = o$observed_TG_neutral

  # get number of observed non-neutral
  elements$observed_non_neutral = 0
  elements$observed_non_neutral[match(o$region_id, elements$region_id)] = o$observed_non_neutral

  # this will be used later to flag elements with many low quality variant calls
  elements$observed_neutral_low_qual = 0
  elements$observed_neutral_low_qual[match(o$region_id, elements$region_id)] = o$observed_neutral_low_qual

  # this will be used later to flag elements with many low quality variant calls
  elements$observed_non_neutral_low_qual = 0
  elements$observed_non_neutral_low_qual[match(o$region_id, elements$region_id)] = o$observed_neutral_low_qual

  # total counts (neutral + non_neutral)
  elements$observed_total = elements$observed_neutral + elements$observed_non_neutral

  # add observed/expected
  #elements$obs_exp_ratio_neutral = elements$observed_neutral/elements$expected_neutral
  #elements$obs_exp_ratio_non_neutral = elements$observed_non_neutral/elements$expected_non_neutral

  elements$observed_low_qual = elements$observed_neutral_low_qual + elements$observed_non_neutral_low_qual

  elements$low_qual_prop = elements$observed_low_qual/(elements$observed_low_qual + elements$observed_neutral + elements$observed_non_neutral + elements$observed_neutral)
  elements$low_qual_prop[is.na(elements$low_qual_prop)] = 0
  elements$filter[elements$low_qual_prop > 0.5 & !(elements$filter == "PASS")] = paste0(elements$filter, ";", "low_quality_0.5_", study$study_name)[elements$low_qual_prop > 0.5 & !(elements$filter == "PASS")]
  elements$filter[elements$low_qual_prop > 0.5 & (elements$filter == "PASS")] = paste0("low_quality_0.5_", study$study_name)  

  # calculate Z score from observed and expected
  var_Z_score = function(observed, expected){
    Xsq_vals = (observed- expected)^2/expected
    excess = ifelse(observed > expected, 1, -1)
    Z = sqrt(Xsq_vals) * excess
  
    # use trimmed z scores to get standard deviation
    Z_trimmed = Z[ (Z > -5) & (Z < 5)]
    Z_trimmed = Z_trimmed[!is.na(Z_trimmed)]
    Z_sd = sd(Z_trimmed)
  
    # divide ALL Z scores by sd from middle set
    Z_normalized = Z/Z_sd
  
    return(Z_normalized)
  }


  # add Z score
  #elements$z_score_neutral = var_Z_score(elements$observed_neutral, elements$expected_neutral)
  #elements$z_score_non_neutral = var_Z_score(elements$observed_non_neutral, elements$expected_non_neutral)

  # rename all of the columns using study name
  colnames(elements)[!(colnames(elements) %in% starting_cols)] = paste0(colnames(elements)[!(colnames(elements) %in% starting_cols)], "_", study$study_name)
}

print(sprintf("Finished annotating elements with variant counts at %s", Sys.time()))
print(head(elements))

if (length(studies) > 1) {
  # now, produce meta obs/exp (non neutral)
  elements$meta_observed_neutral = rowSums(elements[,grepl("observed_neutral", colnames(elements)) & !grepl("low_qual", colnames(elements))])
  #elements$meta_expected_neutral = rowSums(elements[,grepl("expected_neutral", colnames(elements))])
  #elements$meta_obs_exp_ratio_neutral = elements$meta_observed_neutral/elements$meta_expected_neutral
  #elements$meta_z_score_neutral = var_Z_score(elements$meta_observed_neutral, elements$meta_expected_neutral)
  elements$meta_observed_neutral_indels = rowSums(elements[,grepl("observed_neutral_indels", colnames(elements))])
  elements$meta_observed_non_neutral_indels = rowSums(elements[,grepl("observed_non_neutral_indels", colnames(elements))])

  # now, produce meta obs/exp (non neutral)
  elements$meta_observed_non_neutral = rowSums(elements[,grepl("observed_non_neutral", colnames(elements)) & !grepl("low_qual", colnames(elements))])
  #elements$meta_expected_non_neutral = rowSums(elements[,grepl("expected_non_neutral", colnames(elements))])
  #elements$meta_obs_exp_ratio_non_neutral = elements$meta_observed_non_neutral/elements$meta_expected_non_neutral
  #elements$meta_z_score_non_neutral = var_Z_score(elements$meta_observed_non_neutral, elements$meta_expected_non_neutral)

  elements$meta_observed_CA_neutral = rowSums(elements[,grepl("observed_CA_neutral", colnames(elements)) & !grepl("low_qual", colnames(elements))])
  elements$meta_observed_CG_neutral = rowSums(elements[,grepl("observed_CG_neutral", colnames(elements)) & !grepl("low_qual", colnames(elements))])
  elements$meta_observed_CT_neutral = rowSums(elements[,grepl("observed_CT_neutral", colnames(elements)) & !grepl("low_qual", colnames(elements))])
  elements$meta_observed_TA_neutral = rowSums(elements[,grepl("observed_TA_neutral", colnames(elements)) & !grepl("low_qual", colnames(elements))])
  elements$meta_observed_TC_neutral = rowSums(elements[,grepl("observed_TC_neutral", colnames(elements)) & !grepl("low_qual", colnames(elements))])
  elements$meta_observed_TG_neutral = rowSums(elements[,grepl("observed_TG_neutral", colnames(elements)) & !grepl("low_qual", colnames(elements))])
}


# add de novo mutation counts - only SNVs
dnms = read.delim("~/scratch/MUTMODEL/decode_2017_dnms.hg19_liftover.tsv")

dnms = get_region_id_multi_overlap(dnms, elements)
print(head(dnms))

indel_dnms =  subset(dnms, nchar(as.character(ref)) != 1 | nchar(as.character(alt)) != 1)
dnms = subset(dnms, nchar(as.character(ref)) == 1 & nchar(as.character(alt)) ==	1)

#dnms_by_region_id = ddply(dnms, 'region_id', function(df) data.frame(dnm_count = nrow(df), dnm_count_mat = nrow(subset(df, Phase_combined == 'mother')), dnm_count_pat = nrow(subset(df, Phase_combined == 'father'))))
dnms_by_region_id = ddply(dnms, 'region_id', function(df) data.frame(dnm_count = nrow(df)))
indels_by_region_id = ddply(indel_dnms, 'region_id', function(df) data.frame(dnm_count = nrow(df)))

elements$dnm_count_DECODE = 0
elements$dnm_count_DECODE[match(dnms_by_region_id$region_id, elements$region_id)] = dnms_by_region_id$dnm_count

elements$dnm_count_indels_DECODE = 0
elements$dnm_count_indels_DECODE[match(indels_by_region_id$region_id, elements$region_id)] = indels_by_region_id$dnm_count

#elements$dnm_count_mat = 0
#elements$dnm_count_mat[match(dnms_by_region_id$region_id, elements$region_id)] = dnms_by_region_id$dnm_count_mat

#elements$dnm_count_pat = 0
#elements$dnm_count_pat[match(dnms_by_region_id$region_id, elements$region_id)] = dnms_by_region_id$dnm_count_pat

# add DNMs from YUEN 2017 

autism_dnms = read.csv("~/scratch/MUTMODEL/autism_dnms_yuen_2017.csv", stringsAsFactors = F)
colnames(autism_dnms) = c('sample', 'chr', 'start', 'end', 'ref', 'alt', 'platform', 'type', 'validation')
autism_dnms$chr = paste0("chr", autism_dnms$chr)
autism_dnms = get_region_id_multi_overlap(autism_dnms, elements)

indel_dnms =  subset(autism_dnms, nchar(as.character(ref)) != 1 | nchar(as.character(alt)) != 1)
autism_dnms = subset(autism_dnms, nchar(alt) == 1 & nchar(ref) == 1)
dnms_by_region_id = ddply(autism_dnms, 'region_id', function(df) data.frame(dnm_count = nrow(df)))
indels_by_region_id = ddply(indel_dnms, 'region_id', function(df) data.frame(dnm_count = nrow(df)))

elements$dnm_count_YUEN = 0
elements$dnm_count_YUEN[match(dnms_by_region_id$region_id, elements$region_id)] = dnms_by_region_id$dnm_count

elements$dnm_count_indels_YUEN = 0
elements$dnm_count_indels_YUEN[match(indels_by_region_id$region_id, elements$region_id)] = indels_by_region_id$dnm_count

# output annotated elements
write.table(elements, file = args$output, col.names = TRUE, sep = "\t", row.names = FALSE, quote = FALSE)
