# Load the reference sequence as fasta and variants as vcf.
# The vcf is assumed to have haploid organism.
# Output the sequence with ancestral alleles.
# Process the data for one window.

suppressPackageStartupMessages(require(Biostrings))
suppressPackageStartupMessages(require(VariantAnnotation))
suppressPackageStartupMessages(require(parallel))

# Initiate
args <- commandArgs(T)
folder <- args[1]
gff.file <- args[2]
cpu <- args[3]

## TODO: comment after debug
# folder <- "out"
# gff.file <- "chr22.gff"
# cpu <- 12

read_gff <- function(x) {
  d <- read.table(x)
  stringr::str_match(d$V9, "ID=([0-9a-zA_Z]*)")[, 2]
}

make_ancestor_seq <- function(r) {

  # Initiate
  p1 <- sprintf("%s/%s.ref.fa", folder, r)
  p2 <- sprintf("%s/%s.anc.vcf", folder, r)
  
	# Load reference genome
	ref <- readDNAStringSet(p1)
	# Get the start of window in reference
	ws <- as.integer(stringr::str_match(names(ref), ":([0-9]*)-")[, 2])
	# Get the length of window in reference
	wd <- ref@ranges@width
	# Load vcf  
	vcf <- readVcf(p2, "hg19")
	# Get ranges from vcf data
	gr <- rowRanges(vcf)
	# Loop via all ranges and correct DNA reference if required
	for (i in 1:length(gr)) {
	  gt <- geno(vcf)$GT[i]
	  start <- gr[i]@ranges@start
	  pos <- start - ws + 1
	  if (pos < 1 || pos > wd) next
		nuc1 <- "."
		if (gt == 0) nuc1 <- ref(vcf)[i]
		if (gt == 1) nuc1 <- alt(vcf)[[i]]
		nuc2 <- subseq(x = ref, start = pos, width = 1)
		if (nuc2 != nuc1) {
			if (nuc1 == ".") {
				# Change the reference nucleotide
				subseq(x = ref, start = pos, width = 1) <- "-"
			} else {
				# Change the reference nucleotide
				subseq(x = ref, start = pos, width = 1) <- nuc1
			}
		}
	}

	# Set name
	names(ref) <- "ANC"
	# Save ancestor sequence
	p3 <- sprintf("%s/%s.anc.fa", folder, r)
	writeXStringSet(x = ref, filepath = p3)
}

# Get list of region ids
regions <- read_gff(gff.file)
# Loop over all regions in parallel and make ancestor sequences
mclapply(regions, make_ancestor_seq, mc.cores = cpu)



