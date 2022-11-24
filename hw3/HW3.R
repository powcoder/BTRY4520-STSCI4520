#########################################################
#############   BTRY/STSCI 4520  ########################
#############      Homework 3    ########################
############# Due: March 23, 2018 ########################
#########################################################


# Instructions: save this file in the format <NetID>_HW3.R. 
# Complete each question using code below the question number.
# You need only upload this file to CMS. 

# Note, we assume your working directory contains any files
# that accompany this one. 

# Further note: 10% will be deducted if your file produces
# an error when run. If your code produces an error and you 
# cannot find it, comment out that portion of the code and we
# will give partial credit for it. 

# Do not use either the function set.seed() or rm(list=ls())
# in your code. 



#################################
# Question 1: Permutation Tests #
#################################

# The data for this question come from a study on the efficacy
# of medical patch designed to supply a drug to the bloodstream. 
# In this case, the patch was supposed to increase the levels of 
# a hormone. 
#
# The data in 

patchdata = read.table('patchdata.dat',head=TRUE,sep=',')

# contain three measurements for each patient. 
# 
# 1. Placebo measurments
# 2. Old patch -- from an existing manufacturing plant
# 3. New patch -- from a new manufacting plant. 
#
# For this question we will only consider the Placebo
# and old patch measurements, but we will revisit this
# in Question 4. 

# Here we are interested in whether there is any relationship
# between the placebo and patch measurements. To test this,
# we'll make use of the Distance Correlation (see Ch 8 in 
# Rizzo -- Maria Rizzo is one of the inventors of this). 
#
# The idea is to look at pairs of distances between points. 
# If the distance between points in x is similar to the distance
# between the same pair of points in y, then there is a relationship.
#
# Unlike standard correlation, distance correlation will test 
# nonlinear relationships. 

# Formally, we construct the distance correlation as follows. 
#
# First, form the matrix of distances between points
#
#    A_ij =  |x_i - x_j|,    B_ij = |y_i - y_j|
#
# Now obtain the correlation between A and B.  (Technically, 
# we should remove the mean from both row and column, but we'll 
# us a simple version here. 

# a) Write a function to calculate the distance matrix for a 
# vector

distMat = function(x)
{


}

# Using this, write a function for the distance correlation
# between x and y. 

dist.cor = function(x,y)
{

}

# b) Observe that if 

x = seq(-1,1,0.2)
y= x^2

# then cor(x,y) = 0 but dist.cor is not zero, so we are 
# measuring more than the standard correlation. (This is a good
# example to think about why distance correlation works).  

# What is the distance correlation between hormone measurements
# between placebo and old patch in our data?

patch.dist.cor = 

# c) Write a function to conduct a permutation test using 
# distance correlation as a test statistic. You should return
# a list containing elements
#
# - Statistic: the statistic on the original data
# - Null dist: a vector of nperm permuted test statistics
# - p-value: the p-value of the test. 

# Is the test significant for our data? Use nperm = 1000 
# permutations and a significance level of 0.05

perm.dist.cor = function(x,y,nperm)
{


  return(list(obs.correlation = ,
              null.dist = ,
			  pval =  ) )
}


patch.perm = 



# d) We argued that a permutation test always has the
# right level. Let's test this; we'll simulate data
# from the null distribution by making both x and y to
# be standard normals of length 10 (obtained from rnorm(10)). 

# Simulate 1000 data sets and apply the permutation test for
# distance correlation to each using 100 permutations each time.
# We'll use a function wrapper, it should take in the length
# of the vectors, n, the number of permutations nperm and the
# number of simulations nsim and the target level of the 
# test (defaulting to 0.05), return the actual level of the test

Q1dfunction = function(n,nperm,nsim,level=0.05)
{


}



permutation.level = Q1dfunction(10,100,1000)

# e) Repeat the exercise above but let x = rnorm(10) and 
# y = x^2 + rnorm(10). Report the power of the distance 
# correlation test, and the power of the built-in function
# cor.test. 

Q1efunction = function(n,nperm,nsim,level=0.05)
{


  return( c( dist.cor.power, cor.test.power )  )
}



######################################
# Question 2:  False Discovery Rates #
######################################

# In this question we will investigate some
# further properties of false discovery rates

# a) Write a function to simulate a data set
# and carry out an fdr analysis similar to Lecture 
# 8.  It should take inputs
#
#  ngene -- a number of genes
#  nreal -- a number of non-null effects
#  size  -- the difference between means for the 
#           real effects
#  npatients -- number of patients in each group
#  
# and simulate a (2*npatients) x ngene matrix of 
# standard normals, and add size to the first 
# nreal genes of the first npatients.
#
# You should report a vector of p-values and a
# vector of q-values for each gene. These should
# be in the order of the original genes. 
#
# Throughout this question we will set npatients = 10
# and use ngene = 100 and nreal = 20 unless otherwise
# specified. 

fdr.data = function(ngene,nreal,size,npatients)
{


	return( list( p.values = , 
	              q.values = ) )
}


# b) Technically, the false discovery rate is the expected
# percentage of false discoveries (super technically, the
# expectation is less than our estimate), but how variable 
# is it? 

# We will conduct a small simulation study to investigate 
# the sample properties of fdr. For each of the following,
# provide answers based on 100 simulations.

# Write a function to conduct a simulation to generate data 
# from part a nsim times controling the false discovery rate
# at level Q by selecting all genes with q-value less than Q
# and report

# - The percentage of time each gene was selected
# - The number of genes selected each simulation
# - The false discovery proportion each simulation. Set this
#   to zero if your list is of length zero. 
# - The estimated false discovery rate when controlling 
#   pvalues at level P.  (Meaning we expect ngene*P
#   null discoveries.)


fdr.sim = function(ngene,nreal,size,npatients,nsim,Q,P)
{


	return( list( gene.select =  , 
	              num.genes.selected = ,
				  fdp = ,
				  fdr =     ) )
}


# a) What happens if all the genes are null? Use your function
# above to give the mean false discovery proportion and number
# of selected genes when there are no real effects and we control
# at Q = 0.1.  Set ngene = 100, npatients = 10, and P = 0.1.

null.fdp.mean = 

null.no.genes = 

# b) Now we'll allow some real effects. Set nreal=20 and size=3. 
# Here we will be interested in how the total number of tests
# affects the stability of the proportion of false discoveries. 

# Run a simulation with ngene = 100 and with ngene = 1000. Write
# a function that conducts a simulation with nsim data sets and
# reports 
#
#    - the mean and standard deviation of the false discovery
#      proportion (true percentage of wrong discoveries) when the
#      estimated fdr is controlled at Q.  Record these in 
#      mean.fdp and sd.fdp

#    - the mean and standard deviation of your estimated fdr when 
#      you choose all genes with a p-value less than P. Record
#      these in mean.dfr and sd.fdr
#
# 

fdr.sim2 = function(ngene,nreal,size,npatients,nsim,Q,P)
{


	return(list(mean.fdp = 
	            sd.fdp = 
				mean.fdr = 
				sd.fdr =     ) )
}

# Report your results with the following configurations:
#
# 1. ngene = 100, nreal = 20, P = 0.1  this is a baseline
# 2. ngene = 1000, nreal 200, P = 0.1
# 3. ngene = 1000, nreal = 20, P = 0.01
# 4. ngene = 1000, nreal = 20, P = 0.1
#
# For all of these npatient = 10 and Q = 0.1 

# What patterns do you see? Why? Explain how having more
# genes (either real or not) changes the expected number 
# you will pick up and how you expect it to change the 
# variability of the fdr. 





# c) How powerful is this procedure? We'd like to have
# an idea of how likely we are to pick up true results. 
# In our simulation framework, each real difference is
# pretty much the same, so we can treat them as equivalent. 
#
# Write a function to use the output of fdr.sim() and report
#
# 
#  - the percentage of time real effects were detected when fdr
#    is controlled at 0.1 (averaged over the nreal genes).
#  - the percentage of time at least 10 real genes were detected. 
#

Q2cfunction = function( fdr.sim.output, nreal )
{

	return( list( real.one = ,
	              real.ten = ,   ) )
}

# Run a simulation with 20 real genes out of 100 and sets the size 
# parameter to be the values 0.1,0.2,0.4,0.6 taken in turn. Report the two 
# results above in the following vectors:


percent.one = 

percent.ten = 



################################
# Question 3: Multiple Testing #
################################

# This question is associated with recent controversies
# over poor statistical practice, particularly in 
# psychology and nutrition. See a piece from February 25
# in buzzfeed at:
#
# https://www.buzzfeed.com/stephaniemlee/brian-wansink-cornell-p-hacking?utm_term=.nxoJ87k0pG#.ilJpwEGJWX
#
# and a somewhat less aggressive view of the report at
# 
# http://andrewgelman.com/2018/02/27/no-researchers-not-typically-set-prove-specific-hypothesis-study-begins/
# 
# The article describes e-mails that discuss routinely scanning
# many possible hypotheses looking for an effect. 
#
# Here we'll examine one particular case study of an attempt to 
# find significance of a relationship between eating and TV watching 
# in which 400 mediation analyses were trialled. 
#
# In mediation analysis, one seeks to find an effect that is partly 
# obscured by a different covariate. We'll be a bit simpler, and 
# look at trying to find a relationship between y and x while
# controling for another possible covariate z, where there are 400 
# possible z's to look at. 

# We'll consider the following setup. 10 subjects are assigned to 
# watch TV while eating pizza and 10 subjects assigned to read. We
# record the total caloric intake of each subject along with 
# 400 survey variables about their demographics etc. 

# Here y is caloric intake, x is a binary indicator of whether
# the subject watched TV and we have 400 z_i's.  We'll simulate data
# where y has no relationship to x, but might be related to z_i. That
# is we set

x = rep( c(0,1), 1, each = 10)

# and generate y as standard normal 20. We'll generate each z 
# column by
# 
# z = rnorm(1,sd=0.5)*y + rnorm(20)
#
# The first random number here is a 'coefficient' for a relationship 
# with y, except we will end up treating y as a response. 
 
# a) Write a function to generate data as above, run a linear 
# regression y ~ x + z_i  and report the p-value of the 
# coefficient for x for i = 1,...,400. 

Q3a.function = function(){


}


# b) Conduct 100 simulations to determine 
# 
#  i)  How often at least one significant effect is found for x
#  when conducting all tests at the 0.05 level. 
# 
#  ii) How often at least one significant effect is found for x
#  when making a Bonferroni correction.
#
#  iii) The average size of the list of confounders that make 
#  x significant if you apply a false discovery rate procedure 
#  and control the FDR at 0.1.  To do so, obtain the p-values for 
#  x and use this to produce q-values; you can then select the 
#  confounders with q-value less than Q. 
#
# Write your procedure in the following function

Q3simfunc = function(nsim)
{



  return( list( one.significant = ,
                bonferroni.significant = ,
				fdr.size =  ) ) 
}
 
q3result = Q3simfunc(1000)
 
# c) What do you believe is the correct procedure to carry 
# out in this case?



#########################
# Question 4: Bootstrap #
#########################

# In this question we will return to the patch data
# from question 1. A central statistic in evaluating
# the move to a new patch is that the average change
# between the patches should be no larger than 20% of
# the difference between the old patch and the placebo. 
#
# That is we are interested in quantity:
#
# (average[new] - average[old])/(average[old]-average[placebo])
#
# However, ratios are particularly difficult to estimate
# statistically. 

# a) Write a function to obtain a bootstrap distribution for 
# these data. Report the observed statistic, a vector of 
# bootstrap statistics and the 25% and 75% percentiles of the 
# bootstrap distribution. 

ratio.boot = function(data,nboot)
{

 return( list( obs.stat = ,
               boot.stats = , ) )
}


patch.25 = 
patch.75 = 

# b) Write a function to take the result of ratio.boot
# and obtain an estimate of the bootstrap bias. 

ratio.bias = function( boot.obj )
{


}

patch.bias = 


# c) Write a function to take the result of ratio.boot
# and obtain a nonparametric 95% confidence interval for 
# the ratio. 

ratio.confint = function( boot.obj )
{

}

patch.confint = 

# does this confidence interval fall within the required 
# bounds?



##############################
# Question 5: Knockoff Tests #
##############################

# We've seen that you generally can't conduct permutation 
# tests for individual covariates in multiple regression. 
# In one case, however, this is possible. The data in 

peanuts = read.table('peanuts.txt',sep=',',head=TRUE)

# contains data on the yield of peanuts (in pounds) from
# each plant. It's height and whether it had been treated
# with a control,  fast-release fertilizer or slow-release
# fertilizer. 

# Typically, we'd be interested in the effect of treatment,
# but here we will test whether there is a relationship 
# between height and yield, while controling for 
# fertilizer. 

# The key is that we can avoid breaking the relationship
# between fertilizer and height by permuting height 
# WITHIN a level of fertilizer. This will break any 
# relationship between height and yield while keeping 
# the relationship between the covariates. 

# Write a function to conduct a permutation test for 
# the coefficent of height while controlling for 
# fertilizer using the framework above. Report 
# the observed t-statistic, the p-value and the 
# critical value for your test.  We will use
# the estimated coefficient of Height as a statistic.

PartialPermutation = function(data,nperm)
{


   return( list( obs.coef = ,
                 pvalue = ,
				 crit.value = ) )
}

# Test this for yourself using nperm = 1000. 

# This is a special case of the idea of knock-offs
# suggested in 2014 in 
#
#  https://arxiv.org/abs/1404.5609
#
# which were designed for False Discovery Rates. The
# basic idea is that you need to create a new version 
# of x_i that has the same relationship to the other
# covariates, but no relationship to y.  In knockoffs,
# you generate a new  x*_i|x_{-i}   (meaning x_i given 
# all the other covariates). This x*_i has no more 
# information about y than x_{-i} does and so we can
# compare how much it appears to predict to how much
# the real x_i does.  


