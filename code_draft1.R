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
  
# Modelling Decision Tree:
mytree <- rpart(malignant ~ tumor.thickness+uniform.cell.size+uniform.cell.shape+margin.adhesion+
                  epitelial.cell.size+bare.nuclei+bland.chromatin+normal.nucleoli+mitoses+malignant, data = cancer, 
                method = "class") #, control=rpart.control(minsplit=50, cp=0.013))
# Plotting Decision Tree:
rpart.plot(mytree, type = 1, extra=1)
