library(data.table)
library(ggplot2)
library(ggrepel)

samples <- c("501mel", "SKMEL-5", "SKMEL-147")
main_dir <- "ABC"
thres <- 1.2

for (sample in samples) {
  input_path <- file.path(main_dir, sample, "predictions_filtered", "output.tsv")
  output_dir <- file.path(main_dir, sample, "predictions_filtered")
  output <- fread(input_path)
  enhancers <- output[, .(Ngenes = .N, ABC = sum(ABC.Score), abc = mean(ABC.Score)), by = name]
  enhancers[, zscore := resid(lm(ABC ~ Ngenes), type = "pearson")]
  fwrite(enhancers, file.path(output_dir, "enhancer_scores.csv"))
  
  pdf(file.path(output_dir, "ABC_Score_Ranking.pdf"), width = 7, height = 5)
  enhancers[order(ABC), 
            plot(ABC, type = 'b', main = "Enhancer Score Ranking", xlab = "Rank", ylab = "ABC Score", lwd = 2, frame.plot = F)]
  abline(h = thres, col = 2, lty = 2, lwd = 2)
  dev.off()

  top <- enhancers[ABC >= thres]
  fwrite(top, file.path(output_dir, "top_enhancers.csv"))

  p <- ggplot(enhancers, aes(Ngenes, ABC)) +
    geom_point(aes(color = zscore > 1, size = zscore > 1), show.legend = F) +
    geom_smooth(method = "lm") +
    ggrepel::geom_text_repel(aes(label = name), data = enhancers[ABC > 5 & Ngenes < 75]) +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
    scale_size_manual(values = c("FALSE" = .5, "TRUE" = 1)) +
    labs(x = "Number of Genes", y = "Total ABC Score", title = "ABC Score vs Genes") +
    theme_classic()
  ggsave(file.path(output_dir, "ABC_Score_Genes.pdf"), plot = p, width = 7, height = 5)

  bedpe <- output[, .(chrom1 = chr, start1 = start, end1 = end, start2 = TargetGeneTSS, name, score = ABC.Score)]
  bedpe[, `:=`(end2 = start2 + 250, start2 = start2 - 250, chrom2 = chrom1)]
  setcolorder(bedpe, c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score"))
  fwrite(bedpe, file.path(output_dir, "output.bedpe"), sep = "\t", col.names = FALSE)

  cat("Processed sample:", sample, "\n")
}
