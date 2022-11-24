#########################################################
################ BTRY/STSCI 4520 ########################
################## Homework 1 ###########################
############# Due: Feb 9, 2018 #########################
#########################################################


# Instructions: save this file in the format <NetID>_HW1.R. 
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



################################
# Question 1: (from JMR 3.9.1) #
################################

# The function f(x) takes the values
# 
# -x^3     if x < 0
# x^2      if 0 < x < 1
# sqrt(x)  if x > 1
#
# Write a function to compute 

fn = function(x)
{


}

# For the values 

x.values = seq(-2,2,by=0.1)

# create a vector y.values that give the corresponding
# values of f(x). Plot y.values against x.values. 


#################################################
# Question 2:  (from JMR 3.9.2-3.9.4 and 5.7.2) #
#################################################

# For this problem we will consider calculating a 
# sum of powers defining the function
#
# h(x,n) = 1 + x + x^2 + x^3 + ... + x^n
#
# a) Write a function to calculate h(x,n) for any x 
# and n using a for loop.

hn.for = function(x,n)
{


}

# b) Re-write the same function to use a while loop
# instead

hn.while = function(x,n)
{


}

# c) Re-write the function again to use no loops and 
# only employ vector operations

hn.vec = function(x,n)
{


}

# (Note that JMR gives some specific values to check this,
# these will be included in the checking script for the 
# homework. 

# d) As a variation, we are interested in calculating only
# the sum of the even powers
#
# h2(x,n) = 1 + x^2 + x^4 + .... + x^(2*floor(n/2))
#
# where floor(n/2) is the largest whole number smaller than n/2,
# that 2*floor(n/2) is the largest even number smaller than n. 
#
# Modify the code used by any of your functions above to 
# calculate this in 

hn.skip = function(x,n,skip=2)
{


}

# BONUS use the skip argument in hn.skip to sum up 
#
# 1 + x^s + x^(2s) + ... 
#
# for any integer s. 



##################################################
# Question 3: Matrix Multiplication the Slow Way #
##################################################

# Recall that for matrices A and B where A has the same
# number columns as B has rows, we can write the matrix
# product 
#
#  C = AB
#
# by saying that 
#
# C_{ij} =  sum_k  A_{ik} B_{kj}

# Write a function that takes in to matrices A and B and 
# returns C and which only ever sums two numbers together
# or multiplies two numbers together. That is, you are 
# restricted to using only for loops and the mathematical 
# operations + and *. 

matprod = function(A,B)
{



}

# Consider the two matrices 


A = matrix(runif(200*200),200,200)
B = matrix(runif(200*200),200,200)

# How does the timing of 

matprod(A,B)

# Compare to the built-in function

A%*%B 


########################################
# Question 4: Tree Heights (JMR 6.5.3) #
########################################

# a) Create the plot requested in this question. The data in
# treegrowth.txt are available on CMS.  You may do part b first
# and then provide the answer for part a. 

# b) Write a function that extracts a particular
# habitat and produces the corresponding plot

habitatplot = function(num,data)
{



}

# c) Write a function that extracts an individual 
# tree (by its number) and returns the slope in a 
# linear regression of height on age:

ageslope = function(treenum,data)
{

}

###########################################
# Question 5: The Game of Life: JMR 5.7.6 #
###########################################

# You will need to write the function

neighbours = function (A,i,j,n)
{
# A is an n*n 0-1 matrix
# calculate the number of neighbours of A[i,j]


}

# You can find the program to run the game of life
# in life.r; the glider gun initialization is in 
# glidergun.r




###################################
# Bonus: JMR 5.4                  #
###################################

# The function below gives misleading outputs

random.sum <- function(n) {
    # sum of n random numbers
    x[1:n] <- ceiling(10*runif(n))
    cat("x:", x[1:n], "\n")
    return(sum(x))
}

x <- rep(100, 10)

# For example, try

show(random.sum(10))
show(random.sum(5))

# Fix this function so that it does what it says it does. Note in comments below
# which lines you fixed and why they would not work.   You may provide your answers
# with this function commented out. 
