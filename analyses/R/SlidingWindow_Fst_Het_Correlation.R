############# Correlating heterozygosity and Fst ################
#
#
##############################################################



# Load Packages -----------------------------------------------------------

library(dplyr)
library(ggplot2)




# Read in Data ------------------------------------------------------------

swa_het_output_east <- read.delim("Heterozygosity/batch_8_final_filtered_east_HETEROZYGOSITY_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150.txt")
head(swa_het_output_east)

swa_het_output_west <- read.delim("Heterozygosity/batch_8_final_filtered_west_HETEROZYGOSITY_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150.txt")
head(swa_het_output_west)


swa_fst_output <- read.delim("EastvWest/batch_8_SWA_eastwest_globalFst_filtered_output_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150_FILTERED.txt", header = TRUE, sep = "\t")
dim(swa_fst_output)
head(swa_fst_output)



# Join Data Frames --------------------------------------------------------
swa_joined <- full_join(swa_fst_output, swa_het_output_east, by = "position")
head(swa_joined)
colnames(swa_joined) <- c("chromosome","position","Fst.Fct","Mean_boostrap_fst","lower_95_fst","upper_95_fst","pvalue_fst","chromosome.e","He.e","Mean_boostrap.e","lower_95.e","upper_95.e","pvalue.e")
head(swa_joined)
swa_joined_select <- select(swa_joined, c("chromosome","position","Fst.Fct","He.e"))

swa_joined_all <- full_join(swa_joined_select, swa_het_output_west)
head(swa_joined_all)
swa_joined_final <- select(swa_joined_all, c("chromosome","position","Fst.Fct","He.e", "Het"))
colnames(swa_joined_final) <- c("chromosome","position","Fst.Fct","He.e", "He.w")
head(swa_joined_final)


length(which(is.na(swa_joined_final$Fst.Fct)))

## remove NAs
swa_joined_final_filtered <- filter(swa_joined_final, !is.na(Fst.Fct) & !is.na(He.e) & !is.na(He.w))
dim(swa_joined_final)
dim(swa_joined_final_filtered)


# Test Correlation --------------------------------------------------------

east.lm <- lm(swa_joined_final_filtered$Fst.Fct ~ swa_joined_final_filtered$He.e)
summary(east.lm)

#Call:
#  lm(formula = swa_joined_final_filtered$Fst.Fct ~ swa_joined_final_filtered$He.e)

#Residuals:
#  Min       1Q   Median       3Q      Max 
#-0.14185 -0.06513 -0.02672  0.04026  0.78291 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                    0.117300   0.003077  38.118  < 2e-16 ***
#  swa_joined_final_filtered$He.e 0.132397   0.024974   5.301 1.22e-07 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Residual standard error: 0.09283 on 3403 degrees of freedom
#Multiple R-squared:  0.008191,	Adjusted R-squared:  0.0079 
#F-statistic: 28.11 on 1 and 3403 DF,  p-value: 1.222e-07




west.lm <- lm(Fst.Fct ~ He.w, data = swa_joined_final_filtered)
summary(west.lm)

#Call:
#  lm(formula = Fst.Fct ~ He.w, data = swa_joined_final_filtered)

#Residuals:
#  Min       1Q   Median       3Q      Max 
#-0.17858 -0.06101 -0.02675  0.03619  0.79806 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 0.089192   0.003446   25.89   <2e-16 ***
#  He.w        0.345334   0.025239   13.68   <2e-16 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Residual standard error: 0.09075 on 3403 degrees of freedom
#Multiple R-squared:  0.05215,	Adjusted R-squared:  0.05187 
#F-statistic: 187.2 on 1 and 3403 DF,  p-value: < 2.2e-16




# Plot --------------------------------------------------------------------
par(mfrow = c(1,2))
ggplot(swa_joined_final_filtered, aes(x = He.e, y = Fst.Fct)) +
  geom_point() +
  xlab("Heterozygosity, East") +
  ylab("Fst") +
  ylim(0,1) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold")) 
#+
#  annotate("text", x = 0.35, y = 0.75, label = paste("R-squared:", summary(east.lm)$r.squared))

ggplot(swa_joined_final_filtered, aes(x = He.w, y = Fst.Fct)) +
  geom_point() +
  xlab("Heterozygosity, West") +
  ylab("Fst") +
  ylim(0,1) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold")) 
#+
#  annotate("text", x = 0.35, y = 0.75, label = paste("R-squared:", summary(west.lm)$r.squared))








