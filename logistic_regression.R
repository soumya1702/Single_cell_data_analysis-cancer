library(dplyr)

# read in data 
gluta = readRDS('/N/project/ibri_collab/Flak/glu_bst.rds')
# Logistic Regression
meta.data <- gluta[[]] # extract meta data

# first make a table with mice identity as row names and clusters as columns
number_perCluster <- table(meta.data$orig.ident,meta.data$seurat_clusters) 

# drop empty rows
number_perCluster <- number_perCluster[rowSums(number_perCluster[])>0,]

# get proportions
proportion_table <- prop.table(number_perCluster,margin=2)
colSums(proportion_table) # check that column sums equal 1
proportion_table <- as.matrix(proportion_table) # convert table to a matrix

# next create a data frame for the logistic regression model
df <- data.frame(samples=row.names(proportion_table)) # this line creates a row of sample identities
print(df) # check that it works by printing the data frame

# next add a diet column (diet info based on prefix of sample identity)
# logistic regression requires ones and zeros - high-fat diet will be stored as "1"
# For the October data, anything with with an "E" is high-fat diet
# For the other data, anything starting with a "C" is high-fat diet 
df <- df %>% mutate(diet = case_when(
  startsWith(samples,"C") ~ "1", 
  startsWith(samples,"P") ~ "0",
  startsWith(samples,"B") ~ "0",
  startsWith(samples,"October.set1C") ~ "0", 
  startsWith(samples,"October.set2C") ~ "0", 
  startsWith(samples,"October.set1E") ~ "1",
  startsWith(samples,"October.set2E") ~ "1"
))
print(df) # check that there are 2 columns now (sample and diet)

# next add a column to distinguish between male and female samples (again, ones and zeros needed for the regression model)
df <- df %>% mutate(sex = case_when(
  endsWith(samples,"female") ~ "1", 
  endsWith(samples,"_male") ~ "0"
))
print(df) # there should be 3 columns now

# next add a column to distinguish brain regions (bst or poa)
df <- df %>% mutate(region = case_when(
  startsWith(samples,"B") ~ "1", 
  startsWith(samples,"P") ~ "0",
  startsWith(samples,"C-B") ~ "1",
  startsWith(samples,"C-P") ~ "0",
  startsWith(samples,"October.set1C_B") ~ "1", 
  startsWith(samples,"October.set1C_P") ~ "0",
  startsWith(samples,"October.set1E_P") ~ "0",
  startsWith(samples,"October.set2C_P") ~ "0",
  startsWith(samples,"October.set2E_P") ~ "0",
  startsWith(samples,"October.set2E_B") ~ "1",
  startsWith(samples,"October.set2C_B") ~ "1",
  startsWith(samples,"October.set2E_P") ~ "0",
  startsWith(samples,"October.set2C_P") ~ "0"
))
print(df) # check that the df contains 4 columns

# combine df with proportions
mydata <- cbind(proportion_table,region=as.numeric(df$region),diet=as.numeric(df$diet),sex=as.numeric(df$sex))
mydata <- as.data.frame(mydata)
head(mydata) # preview with the head function 

# combine df with proportions
mydata <- cbind(proportion_table,diet=as.numeric(df$diet),sex=as.numeric(df$sex))
mydata <- as.data.frame(mydata)
head(mydata) # preview with the head function 

# GLM - the data is finally ready for the glm function to be ran (generalized linear model)
fit_list = list() # create an empty list to store results

# loop through all clusters *note: change the number in the 'for loop' depending on number of clusters
# you can check the number of clusters with the head(mydata) command above

for(clust in as.character(0:12)){   
  fit_list[[clust]] <- glm(mydata[,clust] ~ diet  + sex  + diet*sex ,data= mydata)
}

# display results
for (i in fit_list){
  print(summary(i))
}

# the best way to save the results is to create a table of useful info and save as a csv file
# initiate lists needed to make a table of p-values
pval_list_diet = list()
pval_list_sex = list()
pval_list_dietsex = list()

# same for effect size 
effectSize_list_diet = list()
effectSize_list_sex = list()
effectSize_list_dietsex = list()

# loop to fill lists - looks confusing but no need to mess with this code, it should run
for(i in fit_list){      
  pvals=summary(i)$coefficients[,4]
  effectSizeDiet=summary(i)$coefficients[,1][2]
  effectSizeSex=summary(i)$coefficients[,1][3]
  effectSizeDietSex=summary(i)$coefficients[,1][4]
  pvalDiet=pvals[2]
  pvalSex=pvals[3]
  pvalDietSex=pvals[4]
  pval_list_diet[[length(pval_list_diet) + 1]] = pvalDiet
  pval_list_sex[[length(pval_list_diet) + 1]] = pvalSex
  pval_list_dietsex[[length(pval_list_diet) + 1]] = pvalDietSex
  effectSize_list_diet[[length(effectSize_list_diet) + 1]] = effectSizeDiet
  effectSize_list_sex[[length(effectSize_list_sex) + 1]] = effectSizeSex
  effectSize_list_dietsex[[length(effectSize_list_dietsex) + 1]] = effectSizeDietSex
}

# create a table with clusters as rows and the different p-vals as columns
pvaldf <- data.frame(unlist(pval_list_diet),unlist(pval_list_sex),
                     unlist(pval_list_dietsex))
names(pvaldf) <- c("Diet","Sex","Diet:Sex")

# be sure to change the number depending on clusters 
row.names(pvaldf) <- paste("Cluster",0:12,sep=' ') # change number here

# save (be sure to name a fitting title)
write.table(pvaldf,'/N/project/ibri_collab/Flak/glu_bst_pvals.csv',row.names = T,col.names = T, sep='\t')

# create a table with cluster as rows and the different effect sizes as columns
effectSizedf <- data.frame(unlist(effectSize_list_diet),unlist(effectSize_list_sex),
                           unlist(effectSize_list_dietsex))
names(effectSizedf) <- c("Diet","Sex","Diet:Sex")

# change number here depending on clusters
row.names(effectSizedf) <- paste("Cluster",0:12,sep=' ') 

# save (be sure to name a fitting title)
write.table(effectSizedf,"/N/project/ibri_collab/Flak/glu_bst_coefficients.csv",row.names = T,col.names = T,sep='\t') 

d2 <- read.csv("/N/project/ibri_collab/Flak/glu_bst_coefficients.csv",header = T, sep='\t')
View(d2)
d1 <- read.csv("/N/project/ibri_collab/Flak/glu_bst_pvals.csv",header = T, sep='\t')
View(d1)