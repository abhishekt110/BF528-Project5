difffile <- read.table("/projectnb/bf528/users/dachshund/project_2/project-2-project-2-dachsund/programmer/cuffdiff_out/gene_exp.diff",
                             header = TRUE)

# sort the above data table so that the smallest q_values are at the top 
sorted_diff <- difffile[order(difffile$q_value),]

#Table with top 10 DE genes, name, FPKM, logfc, pval,qval
top10 <- sorted_diff[1:10,]  
top10DE <- top10[c("gene", "value_1", "value_2", "log2.fold_change.", "p_value", "q_value")]

# produce a histogram of the log2.foldchange column for all genes
logFC_histogram <- hist(difffile$log2.fold_change., 
                        breaks = 60, 
                        main = NULL, 
                        ylim = c(0, 20000),
                        xlim = c(-10, 10),
                        xlab = "Log2 fold change", 
                        col = c("red")
                        )

# look at values where significant == 'YES'
sigDE <- subset(difffile, significant == "yes")

# significant genes histogram
sigDE_hist <- hist(sigDE$log2.fold_change.,
                          breaks = 70, 
                          main = NULL, 
                          ylim = c(0, 500),
                          xlim = c(-10, 10),
                          xlab = "Log2 fold change", 
                          col = c("green")
                          )

# positive and negative logFC
posFC <-subset(sigDE,log2.fold_change. > 0 )
negFC <-subset(sigDE,log2.fold_change. < 0 )
# count genes w/ positive log fold change
dim(posFC)
#count genes w/ negative logFC
dim(negFC)

# write the up and down regulated gene names to separate files
upgenes <- write.csv(posFC$gene, "upregulated_genes.csv")
downgenes <- write.csv(negFC$gene, "downregulated_genes.csv")
