#Running neurocombat on EA values 
#install.packages("devtools")
library(devtools)
#install_github("jfortin1/neuroCombat_Rpackage")
library(neuroCombat)

###FOR RESTING STATE (repeat neurocombat code for EA task MSSD)
#first will combine the regional MSSD together 
SPINS <- read.csv("/projects/tsecara/SPINS_ASD_Project2/data/SPINS/RS_regional_MSSD/SPINS_regional_MSSD_RS.csv")
SPASD <- read.csv("/projects/tsecara/SPINS_ASD_Project2/data/SPASD/RS_regional_MSSD/SPASD_regional_MSSD_RS.csv")

#combine the two dataframes
mssd_val <- rbind(SPASD, SPINS)
write.csv(mssd_val, "/projects/tsecara/SPINS_ASD_Project2/data/combined/regional_MSSD/combined_mssd_mean_RS.csv")

#Loading in the data specific to the neurcombat
add_data <- read.csv("/projects/tsecara/SPINS_ASD_Project2/data/combined/demographics/demo_neurocombat.csv") 
add_data <- add_data[(add_data$record_id %in% mssd_val$record_id), ]
mssd_val <- mssd_val[(mssd_val$record_id %in% add_data$record_id), ] 

###### Removing participant id columns 
rownames(add_data) <- add_data$record_id
add_data <- add_data[,-1]
rownames(mssd_val) <- mssd_val$record_id
mssd_val <- mssd_val[,-1]

#Must transpose so that the participants are in the columns and the features are in the rows  
mssd_val <- t(as.matrix(mssd_val)) #first convert to a matrix, and then transpose 

###Ensure that the validation data is in the right form 
class(add_data$group)
add_data$group <- as.factor(add_data$group)
class(add_data$scanner)
add_data$scanner <- as.factor(add_data$scanner)

#FOR COMBINED SPINS AND SPASD
#mod is a design matrix specifying biological covariates that should be protected - here diagnosis, age, sex, and cog variables
modcombat <- model.matrix(~group +sex + age, data=add_data)

#batch is a vector (length should be equal to the number of columns in the data matrix) that specifies the id for the batch, site, or scanner to correct for
mssd_combat <- neuroCombat(dat=mssd_val, batch=c(add_data$scanner), mod=modcombat)

# transpose the harmonized data matrix
mssd_combat_final <- as.data.frame(t(mssd_combat$dat.combat))

write.csv(mssd_combat_final, "/projects/tsecara/SPINS_ASD_Project2/data/combined/regional_MSSD/neurocombat/combined_mssd_mean_RS_FINAL.csv")