library(biomaRt)
library(EDASeq)

# Read in the Annovar-annotated list of SNPs
data = read.table("chr3_annovar_out.txt", sep = "\t", header=T, na.strings = ".", stringsAsFactors = F, quote='')

# Look up HGNC gene symbols using BiomaRt. 
# Drop SNPs that do not yield gene symbols from the analysis, since they can't be tabulated. 
mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
biomart_data = getBM(c("ensembl_gene_id", "hgnc_symbol", "description"), "ensembl_gene_id", data$Gene.ensGene, mart)
tmp = merge(x = data, y = biomart_data, by.x = "Gene.ensGene", by.y = "ensembl_gene_id")
data_clean = tmp[tmp$hgnc_symbol != "",] # dropping variants with no gene symbol

# Calculate gene length using an EDASeq function to allow normalizing mutational burden by gene length. 
# Here, gene length is calculated as the non-overlapping sum of exons. 
lengths = getGeneLengthAndGCContent(data_clean$Gene.ensGene, org = "hsa")
lengths = lengths[,"length"]

data_clean = unique(cbind(data_clean, lengths)) # get rid of some duplicated rows
colnames(data_clean)[21] = "gene_length"

# save a new copy of the variant data with gene symbols and lengths included
write.table(data_clean, file = "chr3_annotated_all.txt", sep="\t", quote = F, na = ".", row.names = F)
