#BSB513-Homework V: Clustering & Dimension Reduction & Classification

#Question a: Performed the Clustering Analysis

library(GEOquery)
gse13276 <- getGEO("GSE13276", AnnotGPL = TRUE) #using the getGEO function, downloaded the data
exprs_gse13276 <- exprs(gse13276[[1]]) #we defined the expression values on our dataset

#we manipulated the rows and cols. We defined the columns as the samples(genes)
t_gse13276 = t(exprs_gse13276) #samples on the columns in t_gse13276 arg

# We used to pearson correlation in as.dist function.
dist_samples <- as.dist(1-cor(exprs_gse13276, method = "pearson"))

# After, used to hclust function and in dist_samples and centroid method
h_samples <- hclust(dist_samples, method = "centroid")

#we created the dendrogram using the as.dendrogram fucntion
#after, visualize the h_samples dendrogram
plot(as.dendrogram(h_samples))

#Question b: Performed the Principal Component Analysis (PCA)

install.packages("ggplot2")  #for created the plots
library(ggplot2)

pca_result <- prcomp(t_gse13276, center = TRUE, scale = TRUE)
summary(pca_result)

plot(pca_result$x[, 1], pca_result$x[, 2], xlab = "PC1", ylab = "PC2")
text(pca_result$x[, 1], pca_result$x[, 2], rownames(t_gse13276), pos = 3)


#Classification with SVM
install.packages("e1071")
library(e1071)

#tumor sample 0, sur_tissue 1
train_data <- rbind(t_gse13276[1:5,], t_gse13276[6:10,])
label <- c(rep(0, 5), rep(1 ,5))

test_data <- rbind(t_gse13276[11,], t_gse13276[12,])

svm_model <- svm(train_data, label)
print(svm_model)
summary(svm_model)

predictions <- predict(svm_model, train_data)  #test with train data
table(predictions, label) #check accuracy

check <- predict(svm_model, test_data)

cm <- table(check)
cm










