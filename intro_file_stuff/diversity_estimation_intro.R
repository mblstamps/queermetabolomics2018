## Intro to diversity analyses using Amy's+Bryan's models
## And looking at relative abundance modeling

##Loading necessary libraries
library(breakaway)
library(DivNet)
library(corncob)
library(ggplot2)
library(phyloseq)
library(tidyverse)
library(magrittr)

## Loading dataset for examples
data(apples)
data("GlobalPatterns")
data(soil_phylo)

## Subsetting water data from global patterns and glomming to order level
water<-GlobalPatterns %>%
  subset_samples(SampleType %in% c("Freshwater", 
                                   "Freshwater (creek)", 
                                   "Ocean", 
                                   "Sediment (estuary)")) %>%
  tax_glom(taxrank="Order")

## Now for a cursory look at the observed richness
observed_c<-sample_richness(water)
## To visualize outputs as a table
oc_sum<-summary(observed_c)

## Visualize in plot form w/ggplotwrapper
plot(observed_c,water,col='SampleType')

## It appears that AQC has highest richness followed by Sediment
## Let's see if there is confounding with sample depth

depth_frame<-data.frame(obs_r=oc_sum$estimate,
                        depth=sample_sums(water),
                        type=(water %>% sample_data %>% get_variable('SampleType')))

ggplot(depth_frame,aes(x=depth,y=obs_r,col=type))+geom_point()+
  theme_bw()

## And aha, our highest obs_r sample has the highest richness

water_break_alpha<-breakaway(water)
plot(water_break_alpha,water,col='SampleType')

## Aha we have error bars and they swamp differences

## Let's switch and zoom in on one example

trrsed<-subset_samples(water,X.SampleID=="TRRsed1")

## Info: Contains 204 taxa by 7 taxonomic ranks
freq_tab<-trrsed %>% otu_table %>% make_frequency_count_table

## Let's fit the breakaway to just this sample
trr_break<-breakaway(freq_tab)

## What's the model
trr_break$model

## Visualize results
plot(trr_break)

## OKAY IMPORTANT: GENERALLY, YOU RUN BREAKAWAY ON A PHYLOSEQ OBJECT

## Also, we can work with the objects
new_plotframe<-summary(water_break_alpha) %>% 
  add_column('SampleName'= water %>% otu_table %>% sample_names)

## If we wanted a more stable model we can specify
cb_estimate <- water %>%
  chao_bunge

plot(cb_estimate,water,color='SampleType')

## Now that we've observed this let's conduct a hypothesis test
## Hypothesis = all sample types have the same richness
## (at the taxonomic level we're working at, so Orders)
alpha_div_test_res<-betta(summary(water_break_alpha)$estimate,
      summary(water_break_alpha)$error,
      make_design_matrix(water,"SampleType"))

## So let's look at this
## Arg 1: Estimates (from breakaway)
## Arg 2: Errors (from breakaway)
## Arg 3: How to partition data for purpose of hypothesis test

## Now let's look at the result
alpha_div_test_res$table
## These variables are the difference from the mean of the
## first partition from the oother ones

## Now let's learn about how to estimate Shannon diversity
## using DivNet as a network-aware method

## Testing to make sure the package works
divnet_test<-divnet(water,tuning='test')

## DivNet is meant to be run in parallel
water_implementation<-divnet(water,ncores=2)

## Examining the outputs we get
names(water_implementation)

## I am interested in using euclidean metric for beta div
## and maybe Shannon for alpha div but I like the idea
## of estimating richness better

## Visualize this output
plot(water_implementation$shannon,water,col="SampleType")

## Visualize against calculating shannon index of observations
plot(water %>%
  sample_shannon,
  water,ylim=c(0,3.5))

## In the calculation of the metric in the observation a spurious
## trend appears to emerge

## Now do a sample-aware implementation of divnet
sample_water<-divnet(water,X='SampleType',ncores=4)

## Now let's take an look
plot(sample_water$shannon,water,col="SampleType")

## This is very cool, it lets us see the mean and variance changes
## between sample types
## Are there significant differences?

testDiversity(sample_water,h0="shannon")
## The results of this
## we anticipate that the shannon div of freshwater is 3,
## however, all three have significantly lower shannon diversity mean

## Now let's pull out our estimated distance matrix for further testing
bray_curtis_dist<-sample_water$`bray-curtis`

## Notice the distance vectors for samples of the same type will
## be the same because this is to compare ECOSYSTEMS NOT SAMPLES

## We can also get variance estimates (uniquely)

simplifyBeta(sample_water,water,"euclidean","SampleType")

## So we get mean, var, conf intervals

## And we can plot this
simplifyBeta(sample_water, water, "bray-curtis", "SampleType") %>%
  ggplot(aes(x = interaction(Covar1, Covar2), 
             y = beta_est)) +
  geom_point() +
  geom_linerange(aes(ymin = lower, ymax = upper)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") + ylab("Estimates of Bray-Curtis distance")

## Now we can see differences in mean and variances

##LOVE. IT. Gotta try w/phytoplankton data

## Now thinking about corncob -- we want to
## model microbial relative abundance
## so we tackle it with
## and individual taxon regression model
## view our data (soil_phylo)

soil_phylo

## So we have ~8K taxa over 119 samples w/5 covariates

## We take a first look at the OTU table to get oriented
otu_table(soil_phylo)[1:3,1:3]

## These are non-normalized read counts which is great

## Now looking at our five sample variables
sample_data(soil_phylo)[1:3,]

## We have "plant" "day amendment" "amendment" "ID" "Day"

## Now we implement the model
## Restrain analysis to DayAmdmt=11 or 21, so only day 1 and day 2 of amendment 1

soil_for_model<-soil_phylo %>%
  subset_samples(DayAmdmt %in% c('11','21')) %>%
  tax_glom("Phylum")

## To run this function you need
# formula = formula for mean
# phi.formula = formula for overdispersion
# data = phyloseq or data.frame

## So how we run it here
soil_model<-bbdml(formula=OTU.1~1,
                  phi.formula=~1,
                  data=soil_for_model)

## That converged SUPER fast
plot(soil_model)

## Visualization of what's going on
## So what we see is the Relative abundances as observed in each 
## sample with their corresponding confidence intervals
## in absolute abundance space

plot(soil_model,AA=TRUE,color="DayAmdmt")

## Now we model w/covariates
day_model<-bbdml(formula=OTU.1~DayAmdmt,
                 phi.formula=~DayAmdmt,
                 data=soil_for_model)

plot(day_model,color="DayAmdmt")

## Now we can see a difference in the estimated mean and variance
## of the relative abundance of Proteobacteria between
## the day amdmts.

## Is that difference significant statistically?

## FIRST STEP: Likelihood ratio test (is the model that includes our covariate more likely than no)
lrtest(soil_model,day_model)

## Do null model first

## Now for parameter interpretation
summary(day_model)

## Now you get the regression-table-style output for mu and phi
## as long as log likelihood

## Extending the analysis to multiple taxa
## Because we will be doing some network optimization, set seed
set.seed(1)

full_analysis<-differentialTest(formula=~DayAmdmt,
                                phi.formula=~DayAmdmt,
                                formula_null=~1,
                                phi.formula_null=~1,
                                data=soil_for_model,
                                fdr_cutoff=0.05,
                                inits=rbind(rep(0.1,4)))

## Now let's break donw the results
## We get-- p values for tests,
## fdr adjusted p values for tests,
## DA - differentially abundant taxa
## DV - differentially variable taxa