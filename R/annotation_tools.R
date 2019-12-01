library(GenomicRanges)

annotate_mu_background <- function(element, mu_background) {
  e = GRanges(seqnames=Rle(element$chr), ranges = IRanges(start = element$start, end = element$end))
  m = GRanges(seqnames=Rle(mu_background$chr), ranges = IRanges(start = mu_background$start, end = mu_background$end))
  
  nearest_scaff = nearest(e,m)
  filters = mu_background$filter[nearest_scaff]
  
  if (any(filters != 'PASS')) {
    warning('Some of the elements passed are in low-quality regions. Consider filtering with filter == "PASS" before annotating with mu_background.')
  }
  
  mu = mu_background$mu_background[nearest_scaff]
  return(mu)
}

add_mutation_model_annotation = function(elements, mark, mark_name, overlap_type) {
  if (any(!grepl("^chr", mark[,1]))) {
    mark[,1] = paste0("chr", mark[,1])
  }
  
  if (overlap_type == "overlap") {
    if ('pos' %in% colnames(elements)) {
      e = GRanges(seqnames=Rle(elements$chr), ranges = IRanges(start = elements$pos, end = elements$pos))
    } else {
      e = GRanges(seqnames=Rle(elements$chr), ranges = IRanges(start = elements$start, end = elements$end))
    }

    bed = GRanges(seqnames=Rle(mark[,1]), ranges = IRanges(start = mark[,2], end = mark[,3]))
    
    o = findOverlaps(e, bed)
    elements[,mark_name] = FALSE
    elements[unique(queryHits(o)),mark_name] = TRUE  # assumes fourth column of mark is the 'column of interest'

  } else if (overlap_type == "nearest") {
    # assumes mark has chr, pos (not an interval) and need to select nearest position to mark
    if ('pos' %in% colnames(elements)) {
      e = GRanges(seqnames=Rle(elements$chr), ranges = IRanges(start = elements$pos, end = elements$pos))
    } else {
      e = GRanges(seqnames=Rle(elements$chr), ranges = IRanges(start = (elements$end + elements$start)/2, end = (elements$end + elements$start)/2))
    }

    bed = GRanges(seqnames=Rle(mark[,1]), ranges = IRanges(start = mark[,2], end = mark[,2]+1))
    n = nearest(e, bed)
    elements[,mark_name] = mark[n,3]  # assumes third column of mark is the 'column of interest'
  }
  return(elements)
}


get_gene <- function(de_novos, genes){
  # assumes that the gene input has at least three columns with chr, start, end as the first three (all others ignored)
  
  if (any(!grepl("^chr", de_novos$chr))) {
    de_novos$chr = paste0("chr", de_novos$chr)
  }
  if (any(!grepl("^chr", genes$chr))) {
    genes$chr = paste0("chr", genes$chr)
  }
  
  if ("end" %in% colnames(de_novos)){ # region instead of de novos - use the first position of the region to get closest gene!
    dn = GRanges(seqnames=Rle(de_novos$chr), ranges = IRanges(start = de_novos$start, end = de_novos$start + 1))
  } else {
    dn = GRanges(seqnames=Rle(de_novos$chr), ranges = IRanges(start = de_novos$pos, end = de_novos$pos + 1))
  }
  
  if ("end" %in% colnames(genes)) {
    g = GRanges(seqnames=Rle(genes$chr), ranges = IRanges(start = genes$start, end = genes$end))
  } else if ("stop" %in% colnames(elements)) {
    g = GRanges(seqnames=Rle(genes$chr), ranges = IRanges(start = genes$start, end = genes$stop))
  } else {
    stop("Elements did not have 'end' or 'stop' in the column names")
  }
  
  # find overlap between denovos and annotated genes
  hits = findOverlaps(dn, g)
  dn_hits_idx = queryHits(hits) # get index of de novos
  hits_idx = subjectHits(hits) # get index of genes

  de_novos = de_novos[dn_hits_idx,]
  de_novos$gene = genes$gene[hits_idx]
  
  return(de_novos)
}


filter_with_bed <- function(de_novos, bed){
  # assumes that the bed input has at least three columns with chr, start, end as the first three (all others ignored)

  if (any(!grepl("^chr", de_novos$chr))) {
    de_novos$chr = paste0("chr", de_novos$chr)
  }
  if (any(!grepl("^chr", bed[,1]))) {
    bed[ ,1] = paste0("chr", bed[ ,1])
  }
  
  if ("pos" %in% colnames(de_novos)) {
    dn = GRanges(seqnames=Rle(de_novos$chr), ranges = IRanges(start = de_novos$pos, end = de_novos$pos))
  } else {
    dn = GRanges(seqnames=Rle(de_novos$chr), ranges = IRanges(start = de_novos$start, end = de_novos$start))
  }
  gen = GRanges(seqnames=Rle(bed[,1]), ranges = IRanges(start = bed[,2], end = bed[,3]))

  # find overlap between denovos and genomicus annotation
  hits = findOverlaps(dn, gen)
  dn_hits_idx = unique(queryHits(hits)) # get index of de novos with genomicus hit

  de_novos = de_novos[dn_hits_idx, ]

  return(de_novos)
}

exclude_with_bed <- function(de_novos, bed){
  # assumes that the bed input has at least three columns with chr, start, end as the first three (all others ignored)

  if (any(!grepl("^chr", de_novos$chr))) {
    de_novos$chr = paste0("chr", de_novos$chr)
  }
  if (any(!grepl("^chr", bed[,1]))) {
    bed[ ,1] = paste0("chr", bed[ ,1])
  }

  dn = GRanges(seqnames=Rle(de_novos$chr), ranges = IRanges(start = de_novos$pos, end = de_novos$pos))
  gen = GRanges(seqnames=Rle(bed[,1]), ranges = IRanges(start = bed[,2], end = bed[,3]))

  # find overlap between denovos and genomicus annotation
  hits = findOverlaps(dn, gen)
  dn_hits_idx = unique(queryHits(hits)) # get index of de novos with genomicus hit

  de_novos = de_novos[-dn_hits_idx, ]

  return(de_novos)
}



get_region_id <- function(de_novos, elements){
  # assumes that the CNE input has at least three columns with chr, start, end as the first three (all others ignored)

  if (any(!grepl("^chr", de_novos$chr))) {
    de_novos$chr = paste0("chr", de_novos$chr)
  }
  if (any(!grepl("^chr", elements$chr))) {
    elements$chr = paste0("chr", elements$chr)
  }

  if ("end" %in% colnames(de_novos)){ # region instead of de novos - use the first position of the region to get closest gene!
    dn = GRanges(seqnames=Rle(de_novos$chr), ranges = IRanges(start = de_novos$start, end = de_novos$start))
  } else {
    dn = GRanges(seqnames=Rle(de_novos$chr), ranges = IRanges(start = de_novos$pos, end = de_novos$pos))
  }

  if ("end" %in% colnames(elements)) {
    e = GRanges(seqnames=Rle(elements$chr), ranges = IRanges(start = elements$start, end = elements$end))
  } else if ("stop" %in% colnames(elements)) {
    e = GRanges(seqnames=Rle(elements$chr), ranges = IRanges(start = elements$start, end = elements$stop))
  } else {
    stop("Elements did not have 'end' or 'stop' in the column names")
  }

  # find overlap between denovos and annotated elements
  hits = findOverlaps(dn, cne)
  dn_hits_idx = queryHits(hits) # get index of de novos
  CNE_hits_idx = subjectHits(hits) # get index of elements


  return(elements$region_id[CNE_hits_idx])
}

get_region_id_multi_overlap <- function(vars, elements){
  # get the region ID in cases where the elements used might be overlapping (for tiling selection test)

  if (any(!grepl("^chr", vars$chr))) {
    vars$chr = paste0("chr", vars$chr)
  }
  if (any(!grepl("^chr", elements$chr))) {
    elements$chr = paste0("chr", elements$chr)
  }

  if ("end" %in% colnames(vars)){ # region instead of de novos - use the first position of the region to get closest gene!
    v = GRanges(seqnames=Rle(vars$chr), ranges = IRanges(start = vars$start, end = vars$start))
  } else {
    v = GRanges(seqnames=Rle(vars$chr), ranges = IRanges(start = vars$pos, end = vars$pos))
  }

  if ("end" %in% colnames(elements)) {
    e = GRanges(seqnames=Rle(elements$chr), ranges = IRanges(start = elements$start, end = elements$end))
  } else if ("stop" %in% colnames(elements)) {
    e = GRanges(seqnames=Rle(elements$chr), ranges = IRanges(start = elements$start, end = elements$stop))
  } else {
    stop("Elements did not have 'end' or 'stop' in the column names")
  }

  # find overlap between denovos and annotated elements
  hits = findOverlaps(v, e)
  v_hits_idx = queryHits(hits) # get index of de novos
  e_hits_idx = subjectHits(hits) # get index of elements

  vars = vars[v_hits_idx,]  # will repeat rows of v if there are overlapping elements
  vars$region_id = elements$region_id[e_hits_idx]

  return(vars)
}
