#Code applied to individual level data in UK Biobank and HUNT.

#Package requirements

require(data.table)
require(sandwich)
require(lmtest)
require(systemfit)
require(car)
require(survival)
require(nlWaldTest)
require(msm)

#Assumes dataset called "analysis" containing data on IID, FID, phenotypes, covariates and PGS.

#Generate the mean sibship value for measured educational attainment and the educational attainment PGS
analysis$FAM_MEAN <- ave(analysis$Education_std, analysis$FID)
analysis$FAM_MEAN_PGS <- ave(analysis$PGS_std, analysis$FID)

#Specify list of outcome phenotypes
phenotypes <- c("Education", "BMI", "Units", "PackYr", "SBP", "Mortality")

#Generate output dataframe containing observational/PGS association estimates for population and within-sibship models.
output <- data.frame(PHEN=phenotypes, OBS_BETA_POP=NA,  OBS_BETA_WS=NA, OBS_SE_BETA_POP=NA, OBS_SE_BETA_WS=NA, OBS_P_BETA_POP=NA, OBS_P_BETA_WS=NA,
					           PGS_BETA_POP=NA, PGS_BETA_WS=NA, PGS_SE_BETA_POP=NA, PGS_SE_BETA_WS=NA, PGS_P_BETA_POP=NA, PGS_P_BETA_WS=NA)
					  


for (i in 1:6) {

#Specifying first column containing output phenotype
j <- i + 15

# Make a matrix with: [FID PHENOTYPE] [individ - family mean ] [family mean]
merge <- data.table(FID = analysis$FID, OUTCOME = as.numeric(unlist(analysis[,j, with = FALSE])), EXPOSURE = analysis$Education_std, EXPOSURE_PGS = analysis$PGS_std,
 FAM_MEAN = analysis$FAM_MEAN, FAM_MEAN_PGS = analysis$FAM_MEAN_PGS, Sex = analysis$Sex, Survival = analysis$Survival, Dead = analysis$Dead)

#Centre exposure values around the sibship mean
merge2 <- merge[,CENTREDEXPOSURE:=EXPOSURE-FAM_MEAN]
merge2 <- merge[,CENTREDEXPOSURE_PGS:=EXPOSURE_PGS-FAM_MEAN_PGS]



# Run unified regression models for all outcomes asides mortality which requires a different model
if (i < 6) 
{
fit1 <- lm(formula = OUTCOME ~ EXPOSURE, data = merge2)
fit2 <- lm(formula = OUTCOME ~ FAM_MEAN + CENTREDEXPOSURE, data = merge2)
fit3 <- lm(formula = OUTCOME ~ EXPOSURE_GRS, data = merge2)
fit4 <- lm(formula = OUTCOME ~ FAM_MEAN_GRS + CENTREDEXPOSURE_GRS, data = merge2)

#Estimate associations between educational attainment and the 4 outcome phenotypes
if (i > 1) 
{
#SUR_OBS
r1 <- OUTCOME ~ CENTREDEXPOSURE + FAM_MEAN
r2 <- OUTCOME ~ EXPOSURE
merge3 <- merge2[complete.cases(merge2)]
fitsur <- systemfit(list(total = r2, wf = r1), data = merge3)
cluster <- c(merge3$FID, merge3$FID)

vcv_matrixSUR <- vcovCL(fitsur, type = "HC1", cluster)
cov_matrix <- matrix(c(vcv_matrixSUR[4,4], vcv_matrixSUR[2,4], vcv_matrixSUR[2,4], vcv_matrixSUR[2,2]), nrow=2, ncol=2)

}

#Estimate associations between educational attainment PGS and the 5 outcome phenotypes
#SUR_PGS
r1 <- OUTCOME ~ CENTREDEXPOSURE_PGS + FAM_MEAN_PGS
r2 <- OUTCOME ~ EXPOSURE_PGS
merge4 <- merge2[complete.cases(merge2$OUTCOME)]
fitsur_MR <- systemfit(list(total = r2, wf = r1), data = merge4)
cluster <- c(merge4$FID, merge4$FID)

vcv_matrixSUR <- vcovCL(fitsur_MR, type = "HC1", cluster)
cov_matrix <- matrix(c(vcv_matrixSUR[4,4], vcv_matrixSUR[2,4], vcv_matrixSUR[2,4], vcv_matrixSUR[2,2]), nrow=2, ncol=2)

}

#Run specific cox ph models for survival 
if (i == 6) {
surv_object <- Surv(time = merge2$Survival, event = merge2$Dead)
fit1 <- coxph(surv_object ~ Sex + EXPOSURE, data = merge2)
fit2 <- coxph(surv_object ~ Sex + FAM_MEAN + CENTREDEXPOSURE, data = merge2)
fit3 <- coxph(surv_object ~ Sex + EXPOSURE_GRS, data = merge2)
fit4 <- coxph(surv_object ~ Sex + FAM_MEAN_GRS + CENTREDEXPOSURE_GRS, data = merge2)
}

# Extract summary data information

output$OBS_BETA_POP[i] <- fit1$coefficients[2]
output$OBS_BETA_WS[i] <- fit2$coefficients[3]
output$PGS_BETA_POP[i] <- fit3$coefficients[2]
output$PGS_BETA_WS[i] <- fit4$coefficients[3]

# save the variance covariance matrix
vcv_matrix1 <- vcovCL(fit1, cluster=merge2$FID)
vcv_matrix2 <- vcovCL(fit2, cluster=merge2$FID)
vcv_matrix3 <- vcovCL(fit3, cluster=merge2$FID)
vcv_matrix4 <- vcovCL(fit4, cluster=merge2$FID)

#Derive the clustered SEs for the total effect and P-values

test_matrix1 <- coeftest(fit1, vcov.=vcv_matrix1)
test_matrix2 <- coeftest(fit2, vcov.=vcv_matrix2)
test_matrix3 <- coeftest(fit3, vcov.=vcv_matrix3)
test_matrix4 <- coeftest(fit4, vcov.=vcv_matrix4)
	
output$OBS_SE_BETA_POP[i] <- test_matrix1[2,2]
output$OBS_P_BETA_POP[i] <- test_matrix1[2,4]
output$OBS_SE_BETA_WS[i] <- test_matrix2[3,2]
output$OBS_P_BETA_WS[i] <- test_matrix2[3,4]
output$PGS_SE_BETA_POP[i] <- test_matrix3[2,2]
output$PGS_P_BETA_POP[i] <- test_matrix3[2,4]
output$PGS_SE_BETA_WS[i] <- test_matrix4[3,2]
output$PGS_P_BETA_WS[i] <- test_matrix4[3,4]		

}


}
