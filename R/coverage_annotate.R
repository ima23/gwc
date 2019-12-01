library(optparse)
library(plyr)
source("/home/pjs90/software/selection2/R/annotation_tools.R")

option_list <- list(
  make_option("--elements", help = "Pass elements to calculate coverage"),
  make_option("--chromosome", default = 'all', help = "Split to run faster"),
  make_option("--out", help = "output elements with coverage"),
  make_option("--bed_file", default=F, action='store_true', help = "use this option if elements is bed file. assumes first three columns are chr, start, stop"),
  make_option("--rds_file", default=F, action='store_true', help = "use this option if elements is an RDS file - used in BRIDGE regulome project")
)

args <- parse_args(OptionParser(option_list=option_list))

# load in regions file with required columns: chr, start, stop
if (args$bed_file) {
  elements <- read.table(args$elements, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(elements)[1:3] = c("chr", "start", "end")
} else if (args$rds_file) {
  elements <- readRDS(args$elements)
} else {
  elements = read.delim(args$elements)
}

#print(head(elements))
print(nrow(elements))

elements$chr = gsub("^chr", "", elements$chr)  # gnomad coverage files do not have chr in front

if (args$chromosome == 'all') {
  chroms = c(seq(1,22), "X", "Y")
} else {
  chroms = args$chromosome
}

elements = subset(elements, chr %in% chroms)

if (!("region_id" %in% colnames(elements))) {
  if ("stop" %in% colnames(elements)) {
    elements$region_id = paste0(elements$chr, ":", elements$start, "-", elements$stop)
  } else if ("end" %in% colnames(elements)) {
    elements$region_id = paste0(elements$chr, ":", elements$start, "-", elements$end)
  } else {
    stop("No column named 'start' or 'end' in the data provided.")
  } 
  no_region_id = TRUE
} else {
  no_region_id = FALSE
}  

for (i in chroms) {
  i = as.character(i)
  print(sprintf("Working on chromosome %s...", i))

  # use tabix to slice coverage file for this set of elements
  e = subset(elements, chr == i)[,c(1,2,3)]
  e = e[order(e$start),]

  if (any(e$end == e$start)) {
    print('Found lines where end == start... assuming 1-based and subtracting 1 from start position for coverage')
    e$start[e$end == e$start] = e$start[e$end == e$start] - 1
  }


  ### coverage in BRIDGE

  print("Working on BRIDGE coverage")

  element_tmp = tempfile(pattern = "elements", tmpdir = "/tmp")
  output_tmp = tempfile(pattern = "elements", tmpdir = "/tmp")

  write.table(e, file = element_tmp, col.names = FALSE, sep = "\t", row.names = FALSE, quote = FALSE)
  system(sprintf("/home/pjs90/software/tabix-0.2.6/bgzip %s", element_tmp))
  system(sprintf("/home/pjs90/software/tabix-0.2.6/tabix -p 'bed' %s.gz", element_tmp))
  system(sprintf("/home/pjs90/software/tabix-0.2.6/tabix ~/scratch/BRIDGE/coverage/150_20150915_%s.base_median.txt.gz -B %s.gz > %s", i, element_tmp, output_tmp))
  coverage = read.table(output_tmp, fill = T)
  colnames(coverage) = c("chr", "pos", "median_coverage")
  
  #system(sprintf("rm %s.gz", element_tmp))
  #system(sprintf("rm %s", output_tmp))

  #print(head(coverage))
  print(head(elements))

  if (all((elements$end - elements$start) == 0)) { # all of the 'elements' are actually single positions, so merge is much faster
    print('single base elements... annotating with coverage')
    elements = merge(elements, coverage[,c('chr', 'pos', 'median_coverage')], by.x = c('chr', 'start'), by.y = c('chr', 'pos'), all.x = T)
    elements$median_coverage[is.na(elements$median_coverage)] = 0
    elements$median_coverage_BRIDGE = elements$median_coverage
    elements$median_coverage = NULL
  } else {
    print('annotating whole elements with median coverage')
    # give a region id to each line of coverage file using multi region_id
    coverage = get_region_id_multi_overlap(coverage, elements)

    # ddply over the region ids to get coverage per element
    element_coverage = ddply(coverage, "region_id", function(df) data.frame(median_coverage = median(df$median_coverage)))

    print(head(element_coverage))

    elements$median_coverage_BRIDGE = 0
    elements$median_coverage_BRIDGE[match(element_coverage$region_id, elements$region_id)] = element_coverage$median_coverage
  }

  print(nrow(elements))

  ### coverage in gnomAD

  print("Working on gnomAD coverage")

  element_tmp = tempfile(pattern = "elements", tmpdir = "/tmp")
  output_tmp = tempfile(pattern = "elements", tmpdir = "/tmp")

  write.table(e, file = element_tmp, col.names = FALSE, sep = "\t", row.names = FALSE, quote = FALSE)

  system(sprintf("/home/pjs90/software/tabix-0.2.6/bgzip %s", element_tmp))
  system(sprintf("/home/pjs90/software/tabix-0.2.6/tabix -p 'bed' %s.gz", element_tmp))
  system(sprintf("/home/pjs90/software/tabix-0.2.6/tabix ~/scratch/gnomAD/coverage/chr%s.gnomad.median_coverage.txt.gz -B %s.gz > %s", i, element_tmp, output_tmp))
  coverage = read.table(output_tmp, fill = T)
  colnames(coverage) = c("chr", "pos", "median_coverage")  

  print(head(coverage))
  print(head(elements))

  if (all((elements$end - elements$start) == 0)) { # all of the 'elements' are actually single positions, so merge is much faster
    print('single base elements... annotating with coverage')
    elements = merge(elements, coverage[,c('chr', 'pos', 'median_coverage')], by.x = c('chr', 'start'), by.y = c('chr', 'pos'), all.x = T)
    elements$median_coverage[is.na(elements$median_coverage)] = 0
    elements$median_coverage_gnomad = elements$median_coverage
    elements$median_coverage = NULL
  } else {
    # give a region id to each line of coverage file using multi region_id
    coverage = get_region_id_multi_overlap(coverage, elements)

    # ddply over the region ids to get coverage per element
    element_coverage = ddply(coverage, "region_id", function(df) data.frame(median_coverage = median(df$median_coverage)))

    print(head(element_coverage))

    elements$median_coverage_gnomad = 0
    elements$median_coverage_gnomad[match(element_coverage$region_id, elements$region_id)] = element_coverage$median_coverage
  }

  print(head(elements))
}

print(nrow(elements))

elements$chr = paste0("chr", elements$chr)

if (no_region_id) {
  elements$region_id = NULL
}

#if (args$rds_file) {
  #saveRDS(elements, file = args$out)
#} else {
  #write.table(elements, file = args$out, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
#}

write.table(elements, file = args$out, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
