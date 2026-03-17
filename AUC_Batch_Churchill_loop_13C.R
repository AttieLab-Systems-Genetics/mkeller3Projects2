##################################################
#
# Attie RC2 DO Ketones Data
#  batch correction and diet effects
# August 18, 2022
# Gary A. Churchill
##################################################

#####
library(lme4)
library(lmerTest)
library(tidyverse)
library(stringr)
library(lubridate)
library(cowplot)


####
setwd("C:/Users/mkeller3/Desktop")
file <- "C_AUCs_raw.csv"
batch_col <- "run_date_13C"
data_start_col <- 8

# load AUC data for 13C-labeled metabolites
data <- read.csv(file)
data$Batch <- data[, batch_col]

data <- data |>
  mutate(GenLit = factor(GenLit),
         Sex = factor(Sex),
         Diet = factor(Diet),
         Diet.name = factor(Diet.name),
         run_date_13C = mdy(run_date_13C),
         Batch = factor(run_date_13C)) |>
  relocate(Batch, .after = batch_col)

data_start_col <- data_start_col + 1


with(data, table(GenLit, Batch))

anova(lm(AUC_OLEATE_label.2 ~ Sex*Diet+Batch, data = data))


ggplot(data, aes(x = Batch, y = AUC_OLEATE_label.2, color = Diet.name)) + 
  geom_point(position = position_jitter(width = 0.1, height = 0)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))




###############################################################
#
#  batch adjustment functions
#
###############################################################


####
standardise <- function(x) {
  return(
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  )
}

####
outlier_test <- function(x, mult = 2) {
  # Takes a vector of numeric values and returns a logical vector of the same
  # length that indicates whether each value is an outlier

  iqr <- IQR(x, na.rm = TRUE)
  ub <- quantile(x, probs = 0.75, na.rm = TRUE) + mult*iqr
  lb <- quantile(x, probs = 0.25, na.rm = TRUE) - mult*iqr

  return(x > ub | x < lb)
}

####
ab_scale <- function(x, a = 0.01, b = 0.99) {
  # This function rescales a vector of values to a specified range

  MIN <- min(x, na.rm = TRUE)
  RANGE <- max(x, na.rm = TRUE) - MIN

  OUT <- ((x - MIN)/(RANGE)) * (b - a) + a
  return(OUT)
}


################
adjust_for_batch <- function(pheno, batch, sex, diet) {
  # This function uses a fully random effects model to estimate and remove batch
  # effects conditioned on individual, timepoint (batch), sex, and diet.

  ORIG_PHENO <- pheno
  PHENO_NAME <- deparse(substitute(pheno))

  # We will save the original min and max of the data so that we can
  # reset the adjusted values to the same range as the original values.
  ORIGINAL_MIN <- min(pheno, na.rm = TRUE)
  ORIGINAL_MAX <- max(pheno, na.rm = TRUE)


  # Standarised (mean = 0, sd = 1) data for model fitting to reduce risk
  # of convergence errors. The mean and SD used to standardise the data
  # must be saved so that the data can be un-standardised afterwards.
  MEAN <- mean(pheno, na.rm = TRUE)
  SD <- sd(pheno, na.rm = TRUE)
  pheno <- (pheno - MEAN)/SD

  # Convert batch, sex, and diet to characters so that they are treated as
  # categorical variables
  batch <- as.character(batch)
  sex <- as.character(sex)
  diet <- as.character(diet)

  DATA <- data.frame(
    pheno = pheno,
    sex = sex,
    diet = diet,
    batch = batch
  )

  # Fit model to estimate batch effects
  MODEL <- lme4::lmer(
    pheno ~ (1|sex) + (1|diet) + (1|batch),
    control = lme4::lmerControl(
      optimizer = 'bobyqa'
    ),
    data = DATA
  )

  VC <- lme4::VarCorr(MODEL)
  if(sum(VC == 0) > 0 | lme4::isSingular(MODEL)) {
    message('Estimated variance of one of the REF terms is near zero')
  }
  print(VC)

  # Pull estimated batch effects from the model
  REF <- lme4::ranef(MODEL)
  BATCH <- REF$batch
  BATCH <- BATCH[batch, ]

  # Remove estimated batch effects
  PHENO <- pheno - BATCH

  # DATA <- DATA |>
  #   mutate(
  #     adj_transformed = PHENO,
  #     BatchEffect = BATCH
  #   )

  # Un-standardise data
  PHENO <- PHENO * SD + MEAN

  # reset to original range
  PHENO <- ab_scale(PHENO, a = ORIGINAL_MIN, b = ORIGINAL_MAX)

  return(PHENO)
}

###############################################################

# make a separate table of adjusted values

dataAdj <- data

for (i in data_start_col:ncol(data)) {

  dataAdj[, i] <- adjust_for_batch(pheno = data[, i],
                                   batch = data$Batch,
                                   sex = data$Sex,
                                   diet = data$Diet)
}

colnames(dataAdj)[data_start_col:ncol(dataAdj)] <- paste0(colnames(dataAdj)[data_start_col:ncol(dataAdj)], ".Adj")


# option to also make a table containing both unadjusted and adjusted values

combined <- cbind(data, dataAdj[data_start_col:ncol(dataAdj)])


ggplot(combined, aes(x = C_batch, y = AUC_GUANOSINE_label.1.Adj, color = GenLit)) +
  geom_point(position = position_jitter(width = 0.1, height = 0)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))


ggplot(combined, aes(x = C_batch, y = AUC_GUANOSINE_label.1.Adj, color = Diet.name, shape = Sex)) +
  geom_point(position = position_jitter(width = 0.1, height = 0), size = 3) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))


ggplot(combined, aes(x = C_batch, y = AUC_GUANOSINE_label.1, color = Diet.name, shape = Sex)) +
  geom_point(position = position_jitter(width = 0.1, height = 0), size = 3) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))



anova(lm(AUC_N.ACETYLTRYPTOPHAN_label.2 ~ Sex * Diet + Batch, data = combined))
anova(lm(AUC_N.ACETYLTRYPTOPHAN_label.2.Adj ~ Sex * Diet + Batch, data = combined))


Adj <- ggplot(combined, aes(x = Diet.name, y = AUC_N.ACETYLTRYPTOPHAN_label.2.Adj, color = Diet.name)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width=0.2, height=0)) +
  facet_wrap(~Sex) + 
  theme_bw(base_size = 18)

NonAdj <- ggplot(combined, aes(x = Diet.name, y = AUC_N.ACETYLTRYPTOPHAN_label.2, color = Diet.name)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitter(width=0.2, height=0)) +
  facet_wrap(~Sex) +
  theme_bw(base_size = 18)

plot_grid(NonAdj, Adj, labels = "AUTO")

# write file of only adjusted values
write.csv(combined, "C_AUCs_Adj.csv")

# write file of both unadjusted and adjusted values
write.csv(dataAdj, "D_AUCs_raw_Adj.csv")
