#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(facetsSuite)
    library(argparse)
    library(dplyr)
    library(ggplot2)
    library(egg)
    library(purrr)
    library(tibble)
})

args = commandArgs(TRUE)
if (length(args) == 0) {
    message('Run run-facets-wrapper.mini.R --help for list of input arguments.')
    quit()
}

parser = ArgumentParser(description = 'Run FACETS and associated output, input SNP read counts from snp-pileup.')

parser$add_argument('-v', '--verbose', action="store_true", default = TRUE,
                    help = 'Print run info')
parser$add_argument('-f', '--counts-file', required = TRUE,
                    help = 'Merged, gzipped tumor-normal output from snp-pileup')
parser$add_argument('-s', '--sample-id', required = FALSE,
                    help = 'Sample ID, preferrable Tumor_Normal to keep track of the normal used')
parser$add_argument('-D', '--directory', required = TRUE,
                    help = 'Output directory to which all output files are written to')
parser$add_argument('-e', '--everything', dest = 'everything', action = 'store_true',
                    default = FALSE, help = 'Run full suite [default %(default)s]')
parser$add_argument('-g', '--genome', required = FALSE,
                    choices = c('hg18', 'hg19', 'hg38'),
                    default = 'hg19', help = 'Reference genome [default %(default)s]')
parser$add_argument('-c', '--cval', required = FALSE, type = 'integer',
                    default = 50, help = 'Segmentation parameter (cval) [default %(default)s]')
parser$add_argument('-pc', '--purity-cval', required = FALSE, type = 'integer',
                    default = 100, help = 'If two pass, purity segmentation parameter (cval)')
parser$add_argument('-m', '--min-nhet', required = FALSE, type = 'integer',
                    default = 15, help = 'Min. number of heterozygous SNPs required for clustering [default %(default)s]')
parser$add_argument('-pm', '--purity-min-nhet', required = FALSE, type = 'integer',
                    default = 15, help = 'If two pass, purity min. number of heterozygous SNPs (cval) [default %(default)s]')
parser$add_argument('-n', '--snp-window-size', required = FALSE, type = 'integer', 
                    default = 250, help = 'Window size for heterozygous SNPs [default %(default)s]')
parser$add_argument('-nd', '--normal-depth', required = FALSE, type = 'integer',
                    default = 35, help = 'Min. depth in normal to keep SNPs [default %(default)s]')
parser$add_argument('-d', '--dipLogR', required = FALSE, type = 'double',
                    default = NULL, help = 'Manual dipLogR')
parser$add_argument('-S', '--seed', required = FALSE, type = 'integer',
                    default = 100, help = 'Manual seed value [default %(default)s]')
parser$add_argument('-l', '--legacy-output', required = FALSE, type = 'logical',
                    default = FALSE, help = 'create legacy output files (.RData and .cncf.txt) [default %(default)s]')
parser$add_argument('-fl', '--facets-lib-path', required = FALSE,
                    default = '', help = 'path to the facets library. if none provided, uses version available to `library(facets)`')

args = parser$parse_args()

# Helper functions ------------------------------------------------------------------------------------------------

# Write out
write = function(input, output) {
    write.table(input, file = output, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
}

# Print run details
print_run_details = function(outfile,
                             run_type,
                             cval,
                             min_nhet,
                             purity,
                             ploidy,
                             dipLogR,
                             flags = NULL,
                             ...) {
    params = c(...)
    
    run_details = data.frame(
        'sample' = sample_id,
        'run_type' = run_type,
        'purity' = signif(purity, 2),
        'ploidy' = signif(ploidy, 2),
        'dipLogR' = signif(dipLogR, 2),
        'facets_version' = as.character(packageVersion('facets')),
        'cval' = cval,
        'snp_nbhd' = args$snp_window_size,
        'min_nhet' = min_nhet,
        'ndepth' = args$normal_depth,
        'genome' = args$genome,
        'seed' = args$seed,
        'flags' = flags,
        'input_file' = basename(args$counts_file))
    
    if (length(params) > 0) {
        run_details = data.frame(run_details,
                                 'genome_doubled' = params$genome_doubled,
                                 'fraction_cna' = signif(as.numeric(params$fraction_cna), 2),
                                 'hypoploid' = params$hypoploid,
                                 'fraction_loh' = signif(as.numeric(params$fraction_loh), 2),
                                 'lst' = params$lst,
                                 'ntai' = params$ntelomeric_ai,
                                 'hrd_loh' = params$hrd_loh)
    }
    
    write(run_details, outfile)
    
    run_details
}

# Default set of output plots
print_plots = function(outfile,
                       facets_output,
                       cval) {
    
    plot_title = paste0(sample_id,
                        ' | cval=', cval,
                        ' | purity=', round(facets_output$purity, 2),
                        ' | ploidy=', round(facets_output$ploidy, 2),
                        ' | dipLogR=', round(facets_output$dipLogR, 2))
    
    png(file = outfile, width = 850, height = 999, units = 'px', type = 'cairo-png', res = 96)
    suppressWarnings(
        egg::ggarrange(
            plots = list(
                cnlr_plot(facets_output),
                valor_plot(facets_output),
                icn_plot(facets_output, method = 'em'),
                cf_plot(facets_output, method = 'em'),
                icn_plot(facets_output, method = 'cncf'),
                cf_plot(facets_output, method = 'cncf')
            ),
            ncol = 1,
            nrow = 6,
            heights = c(1, 1, 1, .15, 1, .15),
            top = plot_title)
    )
    dev.off()
}

# Print segmentation
print_segments = function(outfile,
                          facets_output) {
    write(facets_output$segs, outfile)
}

# Print IGV-style .seg file
print_igv = function(outfile,
                     facets_output) {
    
    ii = format_igv_seg(facets_output = facets_output,
                        sample_id = sample_id,
                        normalize = T)
    
    write(ii, outfile)
}

# Define facets iteration
# Given a set of parameters, do:
# 1. Run facets
# 2. Generate and save plots
# 3. Print run iformation, IGV-style seg file, segmentation data
facets_iteration = function(name_prefix, ...) {
    params = list(...)
    
    output = run_facets(read_counts = read_counts,
                        cval = params$cval,
                        dipLogR = params$dipLogR,
                        ndepth = params$ndepth,
                        snp_nbhd = params$snp_nbhd,
                        min_nhet = params$min_nhet,
                        genome = params$genome,
                        seed = params$seed,
                        facets_lib_path = params$facets_lib_path)
    
    
    output
}

# Run -------------------------------------------------------------------------------------------------------------

# Name files and create output directory
sample_id = ifelse(is.na(args$sample_id),
                   gsub('(.dat.gz$|.gz$)', '', basename(args$counts_file)),
                   args$sample_id)
directory = args$directory

if (dir.exists(directory)) {
    #stop('Output directory already exists, specify a different one.',  call. = F)
} else {
    system(paste('mkdir -p', directory))
}

# Read SNP counts file
message(paste('Reading', args$counts_file))
read_counts = read_snp_matrix(args$counts_file)
message(paste('Writing to', directory))

# Determine if running two-pass
if (FALSE) {
    
} else {
    name = paste0(directory, '/', sample_id)

    output = facets_iteration(name_prefix = name,
                              dipLogR = args$dipLogR,
                              cval = args$cval,
                              ndepth = args$normal_depth,
                              snp_nbhd = args$snp_window_size,
                              min_nhet = args$min_nhet,
                              genome = args$genome,
                              seed = args$seed,
                              facets_lib_path = args$facets_lib_path)
    saveRDS(output, paste0(directory, '/', sample_id, '.rds'))
}
