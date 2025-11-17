library(GENOVA)
library(rhdf5)
library(strawr)

setwd('Plots')
dir <- 'Plots'

# Load HiC matrix
SKMEL5_10kb <- load_contacts(signal_path = 'SKMEL-5.mcool',
                                    sample_name = "SKMEL-5",
                                     resolution = 10000,
                                     balancing = T, colour = "red")

MEL501_10kb <- load_contacts(signal_path = '501mel.mcool',
                                     sample_name = "501mel",
                                     resolution = 10000,
                                     balancing = T, colour = "red")

SKMEL5_loops = read.delim('/gpfs/data/HernandoLab/home/bericp01/HiC/results/loops/peakachu/SKMEL-5/SKMEL-5_validate_CTCF_merged_loops_matched_locations', h = F)
SKMEL5_tads = read.delim('/gpfs/data/HernandoLab/home/bericp01/HiC/results/TADs/domaincaller/SKMEL-5/SKMEL-5_40kb.output', h = F)

MEL501_loops1 = read.delim('/gpfs/data/HernandoLab/home/bericp01/HiC/results/loops/peakachu/501mel_siNTC/501mel_siNTC_validate_CTCF_merged_loops_matched_locations', h = F)
MEL501_tads1 = read.delim('/gpfs/data/HernandoLab/home/bericp01/HiC/results/TADs/domaincaller/501mel_siNTC/501mel_siNTC_40kb.output', h = F)

####Region HOXD13
pdf(file.path(dir, "HiC_Genova_pyramid_SKMEL-5_HOXD13.pdf"))
p <- pyramid(SKMEL5_10kb,
  chrom = 'chr2',
  start = 175.0e6,
  end = 177.0e6,
  colour = c(0,30), edge = "black")
p + add_tads(SKMEL5_tads, colour = "red") + add_loops(SKMEL5_loops, colour = "red", shape = 1)
dev.off()

pdf(file.path(dir, "HiC_Genova_pyramid_501mel_HOXD13.pdf"))
p <- pyramid(MEL501_10kb,
  chrom = 'chr2',
  start = 175.0e6,
  end = 177.0e6,
  colour = c(0,30), edge = "black")
p + add_tads(MEL501_tads, colour = "red") + add_loops(MEL501_loops, colour = "red", shape = 1)
dev.off()

pdf(file.path(dir, "HiC_Genova_pyramid_SKMEL5_501mel_diff_HOXD13.pdf"))
pyramid_difference(SKMEL5_10kb, MEL501_10kb,
  chrom = 'chr2', start = 175.0e6, end = 177.0e6)
dev.off()

