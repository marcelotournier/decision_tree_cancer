# loading cancer biopsy data:
library(MASS)
library(rpart)
library(rpart.plot)
# loading all columns, except ID
cancer <- MASS::biopsy[,2:11]
# renamiing cols:
colnames(cancer) <- c("tumor.thickness",
                     "uniform.cell.size",
                     "uniform.cell.shape",
                     "margin.adhesion",
                     "epitelial.cell.size",
                     "bare.nuclei", # 16 NAs here -- we will replace them by medians
                     "bland.chromatin",
                     "normal.nucleoli",
                     "mitoses",
                     "malignant"
                     )
# recoding target(dependent) variable malignant to 0 or 1
cancer$malignant <- as.numeric(as.factor(cancer$malignant))-1
# NA Fills - Checking NA Values:
#View(cancer[is.na(cancer$bare.nuclei),])
# only 2 of 14 values are malignant.  We will fill with the bare.nuclei median value for benign tumors
#filling missing values using vector operations (more efficient than for loops!)
cancer[is.na(cancer$bare.nuclei),]$bare.nuclei <- median(cancer[cancer$malignant == 0,]$bare.nuclei, na.rm = T)

#building train-test splits:
#generating test and train data - Data selected randomly with a 80/20 split
x <- cancer
trainIndex  <- sample(1:nrow(x), 0.8 * nrow(x))
train <- x[trainIndex,]
test <- x[-trainIndex,]
 
# Modelling Decision Tree:
mytree <- rpart(malignant ~ tumor.thickness+uniform.cell.size+uniform.cell.shape+margin.adhesion+
                  epitelial.cell.size+bare.nuclei+bland.chromatin+normal.nucleoli+mitoses+malignant, data = train, 
                method = "class") #, control=rpart.control(minsplit=50, cp=0.013))

# testing predictions on the model:
t_pred <- predict(mytree,test,type="class")
t <- test['malignant']

#building a confusion matrix
confMat <- table(test$malignant,t_pred)

# getting model accuracy:
accuracy <- sum(diag(confMat))/sum(confMat)
print(accuracy) # 0.9285714 -> 92.8%  Not bad!!!

# Plotting our Decision Tree:
rpart.plot(mytree, type = 1, extra=1)


