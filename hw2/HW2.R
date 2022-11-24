#########################################################
#############   BTRY/STSCI 4520  ########################
#############      Homework 2    ########################
############# Due: March 2, 2018 ########################
#########################################################


# Instructions: save this file in the format <NetID>_HW2.R. 
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



##########################################
# Question 1: (from JMR 4.6.1 and 4.6.2) #
##########################################

# a)  4.6.1, using data files from CMS.   
# You should write your code in a function
# that reads in two filenames (as character strings) 
# and outputs a table; make sure the subject id 
# (which you can assume is the first column) matches
# between the data set, and produce an NA when the 
# a subject appears in one data set and not in the
# other. 

merge.files = function(filename1,filename2)
{



}


# b) From 4.6.2; write a function to order a data frame
# by one of its columns. That is, to re-arrange the 
# rows so that a column (indicated by col.ind in the
# input below) are given from smallest to largest. 
# You may use the sort or order functions in R. 
# You can check your result using the teeth data.

order1 = function(data, col.index)
{


}


# c) You will note that if we order the tooth data 
# by age, there are ties at some ages. Write a function
# that will order a function by the entries in column
# col.index1 but will break ties in those entries using 
# column col.index2 and return the resulting sorted 
# data set. 

order2 = function(data, col.index1, col.index2)
{



}

# Use this to sort your data set first by age and
# then by teeth. 



#############################
# Question 2:  An Easy Sort #
#############################

# a) Write a function that takes in two sorted vectors
# a and b, and returns c: a sorted vector of their
# combined entries. 
#
# Do so as efficiently as possible, WITHOUT using
# built-in functions; only element-wise comparisons
# and assignments. 

two.sort = function(a,b)
{

  return(c)
}

# You can check the working of your function using

a = sort(rnorm(1000))
b = sort(rnorm(500))

c = two.sort(a,b)

plot(diff(c))   # These should all be positive. 

# b) If a has N1 elements and b has N2 elements, under
# what configuration do you do the smallest number 
# of comparisons (and how many comparisons do you 
# make)? Under what configuration do you do the 
# largest?
#
# Give your answer in comments. 



# BONUS:  Assuming that the entries of a and b are
# from the same distribution, given an average-case
# order of computational cost in terms of 
# N1 and N2. 





################################################
# Question 3: A Numerical Stability Experiment #
################################################

# This question examines three methods for evaluating
# the polynomial 
#
# f(x) = a_n x^n + a_{n-1} x^{n-1} + ... + a_1 x + a_0
#
# The first uses the expression as given, we'll call this
# 'standard'.
#
# The second of these is Horner's method that represents
# the polynomial as a set of multiplications
#
# f(x) = (...((a_n x + a_{n-1}) x + a_{n-2}) ... a_1) x + a_0
#
# The third is the factorization method if we know the 
# roots of the polynomial are r_1,...,r_n:
#
# f(x) = r_0 (x - r_1)(x-r_2)...(x-r_n)

# a) Write functions standard, horner and factored that evaluate 
# f(x) with arguments x and vectors a or r as appropriate. 
# these functions should not assume a fixed degree n of the 
# polynomial. 

standard = function(x,a)
{

}


horner = function(x,a)
{


}

factored = function(x,r)
{


}

# b) Use these functions to evaluate 
# 
# f(x) = (x-2)^9 
#      = x^9 - 18x^8 + 144x^7 - 672x^6 + 2016x^5
#          - 4032x^4 + 5376x^3 -4608x^2 + 2304x - 512
#
# on a grid of 101 values of x from x = 1.95 to x = 2.05.
#
# Record the results of each evaluation in vectors named
# y.standard, y.horner and y.factored respectively. 
#
# Plot the results of these evaluations in different 
# colors. Which do you feel is most accurate?

y.standard = 

y.horner = 

y.factored = 

# c) Now repeat your plot but for the region x = 1.5 to 2.5;
# do the differences appear to be important? Give a reason for
# the differences between the three methods of evaluating the
# same function.   
#
# Do NOT overwrite the values in y.standard, y.horner and y.factored. 


# d) One consequence of these differences is in finding the roots
# of a function; ie where f(x) = 0 (in this case we know that's at
# x = 2).  We can do this for looking at the places were f(x) < 0
# on one side and f(x) > 0 on the other.  By default, we'll record
# the value on the left side. So if I have x[14] < 0 < x[15], we'll
# record x[14] as the zero. 
#
# Using the evaluations from part b, report the zeros found by 
# the three methods.  Record these in vectors zero.standard,
# zero.horner and zero.factored respectively. 

zero.standard = 

zero.horner = 

zero.factored = 



# BONUS 2 points for constructing your functions to accept a 
# vector of x values and using this functionality throughout. 



###################################
# Question 4: Simulated Censoring #
###################################

# Clinical trials face a problem when the primary endpoint
# is how long patients live after treatment -- there isn't
# the time to wait for everyone to die! This means that at
# the end of the trial (say 5 years), some observations are
# "censored" -- we know that the patient lived at least 5
# years but don't know how much longer. 
#
# Generally, these data are given as (T_i,S_i) where 
# T_i is the observed time of age or drop-out and S_i 
# tells you whether the observation was censored. 
#
# There are some sophisticated ways to deal with this, but
# a naive way to use these data to estimate life expectancy are:
#
# 1. Keep all the data together and just used the 
# censor time (ie 5 years in this case) for the 
# censored observations.
#
# 2. Only look at the non-censored data; ie, those
# who are already dead. 
#
# We know that both of these will under-estimate life
# expectancy but we would like to know by how much. 

# To investigate this, we'll consider post-treatment 
# lifetimes as having an exponential distribution with
# rate 0.3.  You can generate n samples from this 
# distribution with the function rexp(n,0.3). 

# a) Write a function to generate n samples from this
# distribution and return three data sets:
#
# 1. The original, uncensored lifetimes. 
#
# 2. The lifetimes where any that are over 5 years
# are replaced with 5. 
#
# 3. The lifetimes with anything over 5 years removed. 

lifedata = function(n)
{


}

# b) Using your function lifedata, simulate 100 
# data points and calculate the differences between
# the average uncensored lifetimes (1) and the 
# truncated lifetimes (2).  Give an estimate 
# of the standard deviation of this difference. 
#
# Record this mean and standard deviation in objects
# meandiff and sd.meandiff  


meandiff = 
sd.meandiff = 


# c) Calculate how many data points you would need to
# ensure that your estimated difference has standard
# deviation less that 0.01.   Store this in npts.diff

npts.diff = 


# d) Run 100 simulations using the data size you 
# calculated in part c. Is your standard deviation less
# than 0.1?  Calculate the standard deviation of the
# difference between the mean true lifetimes and the
# mean of the lifetimes that are less than 5. 
#
# Record your standard deviations in sd.trunc for the
# difference between original and truncated data, and 
# sd.remove for the difference between the original data
# and when you've removed the censored observations. 

sd.trunc = 
sd.remove = 


# e) Why is it difficult to repeat c) for the difference
# between the mean of the true lifetimes (1) and the
# reduced data set where you only keep those less than 5 (3)?



############################
# Question 5: Gini Indices #
############################

# The Gini ratio is a measure of inequality in a data
# set of n positive numbers (x_1,...,x_n) and is defined 
# by
#
# G = sum_i sum_j | x_i - x_j | / (2 n^2 bar(x) )
#
# Where bar(x) is the average of the x's. 

# a) Write a function that produces the Gini index for
# a given data set using as few loops (including apply 
# statements) as possible. (Bonus if you can use none). 

Gini = function(x)
{


}

# It may help to first obtain the matrix D_ij = |x_i - x_j|. 

# b) Simulate 1000 data sets of size 20 from an Exponential(1) 
# distribution and use these to calculate the mean and variance 
# the Gini index for data generated this way. 
#
# Store these in objects gmean and gvar respectively. 

gmean = 

gvar = 


# c) Consider testing whether the Gini index is larger than 
# Exp(1) data would give you by rejecting if G > mean + 2 sd
# where mean and sd are calculated from your result in part b. 
# When the data are Exp(1), what percentage of the time do you 
# reject the test?  Record this in the obeject gini.alpha.
#
# You should not need to generate new data. 

gini.alpha = 

# d) Write a function to calculate the power of this test for
# any exponential rate parameter mu; ie where you simulate from
# rexp(n,mu).  It should have arguments given by your mean (argument m) 
# and sd (argument s), the number of simulations to run, and mu:

gini.power = function(mu, m, sd, nsim)
{


}


# e) Construct a power curve using your function in part d)
# 20 values of mu from mu = 0.5 to mu = 2.5. Plot this curve. 


# (The following will not be graded)
# You may want to try using an alternative log normal distribution
# with parameter m; you can generate n samples from this with
#
# exp( m*rnorm( n ) )
#
# See what happens when you generate a power curve with this family
# of distributions, taking m to be the values of mu you used above. 


# BONUS: Repeat part b without using for loops. It may help
# to think about how to form the vectorized version of D_ij 
# using matrix multiplication; to do that, you may want to 
# consider kronecker products obtained from %x% in R.  