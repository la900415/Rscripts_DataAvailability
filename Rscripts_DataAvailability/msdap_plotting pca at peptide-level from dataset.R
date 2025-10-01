############ plot pca of peptides from dataset ################################
# peptide level data matrix, using filtered+normalized peptides across all groups
tibw_noexclude = dataset$peptides %>% 
  select(peptide_id, sample_id, intensity_all_group) %>%
  filter(!is.na(intensity_all_group)) %>%
  pivot_wider(id_cols = peptide_id, names_from = sample_id, values_from = intensity_all_group)
matrix_sample_intensities = msdap:::as_matrix_except_first_column(tibw_noexclude)

# compute PCA
PPCA = pcaMethods::pca(base::scale(t(matrix_sample_intensities), center=TRUE, scale=FALSE), method="ppca", nPcs = 3, seed = 123, maxIterations = 2000)
mat_pca = PPCA@scores # PCA coordinate matrix; rownames = sample_id and colnames = PCA dimension
pca_var = array(PPCA@R2, dimnames = list(colnames(PPCA@scores))) # variation per dimension

# if you only want a table with the PCA coordinates, just print the mat_pca variable and you're done at this point

# example plot code
dims = c(1,2) # plot dimensions. in this example, can be; 1:2, 2:3, c(1,3)
p = ggplot(data.frame(x=mat_pca[,dims[1]], y=mat_pca[,dims[2]], label=rownames(mat_pca)), 
           aes(x = x, y = y, label = label)) +
  geom_point() +
  geom_text() +
  labs(x = sprintf("dimension %d (%.1f%%)", dims[1], pca_var[dims[1]] * 100),
       y = sprintf("dimension %d (%.1f%%)", dims[2], pca_var[dims[2]] * 100) ) +
  theme_bw()
print(p)