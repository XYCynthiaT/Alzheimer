setwd(dirname(parent.frame(2)$ofile))

mayo <- readRDS("../data/mayo.rds")
msbb <- readRDS("../data/msbb.rds")
rosmap <- readRDS("../data/rosmap.rds")

library(tidyverse)


table(mayo$TCX$pdata$Diagnosis, mayo$TCX$pdata$Gender)
table(mayo$TCX$pdata$Diagnosis, mayo$TCX$pdata$AgeAtDeath)
ifelse(is.na(mayo$TCX$pdata$ApoE), "NA", 
                ifelse(mayo$TCX$pdata$ApoE %in% c(24, 34, 44), "present", "absent")) %>%
    table()
ifelse(is.na(mayo$TCX$pdata$ApoE), "NA", 
       ifelse(mayo$TCX$pdata$ApoE %in% c(24, 34, 44), "present", "absent")) %>%
    table(mayo$TCX$pdata$Diagnosis)
table(mayo$TCX$pdata$Diagnosis, mayo$TCX$pdata$ApoE)
table(mayo$TCX$pdata$Diagnosis, mayo$TCX$pdata$race)
table(mayo$TCX$pdata$Diagnosis)
table(mayo$TCX$pdata$Gender)
table(mayo$TCX$pdata$race)

table(rosmap$dorsolateral$pdata$diagnosis, rosmap$dorsolateral$pdata$msex)
table(rosmap$dorsolateral$pdata$diagnosis, rosmap$dorsolateral$pdata$race)
ifelse(is.na(rosmap$dorsolateral$pdata$apoe_genotype), "NA", 
       ifelse(rosmap$dorsolateral$pdata$apoe_genotype %in% c(24, 34, 44), "present", "absent")) %>%
    table()
ifelse(is.na(rosmap$dorsolateral$pdata$apoe_genotype), "NA", 
       ifelse(rosmap$dorsolateral$pdata$apoe_genotype %in% c(24, 34, 44), "present", "absent")) %>%
    table(rosmap$dorsolateral$pdata$diagnosis)
table(rosmap$dorsolateral$pdata$diagnosis)
table(rosmap$dorsolateral$pdata$msex)
table(rosmap$dorsolateral$pdata$race)


for (area in unique(msbb$pdata$BrodmannArea)) {
    data <- subset_samples(msbb, msbb$pdata$BrodmannArea==area & msbb$pdata$diagnosis!="Transition") 
    print(area)
    print(table(data$pdata$diagnosis, data$pdata$sex))
}
for (area in unique(msbb$pdata$BrodmannArea)) {
    data <- subset_samples(msbb, msbb$pdata$BrodmannArea==area & msbb$pdata$diagnosis!="Transition") 
    print(area)
    print(table(data$pdata$diagnosis, data$pdata$race))
}
for (area in unique(msbb$pdata$BrodmannArea)) {
    data <- subset_samples(msbb, msbb$pdata$BrodmannArea==area & msbb$pdata$diagnosis!="Transition") 
    print(area)
    print(table(data$pdata$diagnosis))
}
for (area in unique(msbb$pdata$BrodmannArea)) {
    data <- subset_samples(msbb, msbb$pdata$BrodmannArea==area & msbb$pdata$diagnosis!="Transition") 
    apoe4 <- ifelse(is.na(data$pdata$apoeGenotype), "NA", 
                    ifelse(data$pdata$apoeGenotype %in% c(24, 34, 44), "present", "absent"))
    print(area)
    print(table(data$pdata$diagnosis, apoe4))
}
msbb_sub <- subset_samples(msbb, msbb$pdata$diagnosis!="Transition") 
table(msbb_sub$pdata$sex)
table(msbb_sub$pdata$race)
apoe4 <- ifelse(is.na(msbb_sub$pdata$apoeGenotype), "NA", 
                ifelse(msbb_sub$pdata$apoeGenotype %in% c(24, 34, 44), "present", "absent"))
table(apoe4)