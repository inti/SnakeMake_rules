library(readr)
library(tximport)

tx2gene <- read_csv("~/DATA/ref_genomes/B_terrestris/Bombus_terrestris.Bter_1.0.cdna.all.fa.gene_transcript.txt")
samples = read_table("SRR_Acc_list.txt",col_names="run")
files <- file.path("rnaseq_quant/salmon", samples[["run"]], "quant.sf")
names(files) = samples[["run"]]
tx_quant <- tximport(files, type = "salmon", tx2gene = tx2gene, txOut=TRUE)
gene_quant <- summarizeToGene(tx_quant, tx2gene,countsFromAbundance="lengthScaledTPM")

write.table(as.data.frame(gene_quant$counts), file = snakemake@output[["gene_lengthScaledTPM"]], sep=" ",quote=F)
write.table(as.data.frame(gene_quant$counts), file = snakemake@output[["gene_length"]], sep=" ",quote=F)

gene_quant <- summarizeToGene(tx_quant, tx2gene)
write.table(as.data.frame(gene_quant$counts), file = snakemake@output[["gene_counts"]], sep=" ",quote=F)

write.table(as.data.frame(tx_quant$counts), file = snakemake@output[["tx_counts"]], sep=" ",quote=F)
write.table(as.data.frame(tx_quant$counts), file = snakemake@output[["tx_length"]], sep=" ",quote=F)

print(snakemake@output)
