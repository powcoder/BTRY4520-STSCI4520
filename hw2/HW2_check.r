https://powcoder.com
代写代考加微信 powcoder
Assignment Project Exam Help
Add WeChat powcoder
https://powcoder.com
代写代考加微信 powcoder
Assignment Project Exam Help
Add WeChat powcoder
##############################
##### Homework 2 CheckFile####
##############################

# First we'll call your solutions
# if these are giving you an error, you
# can add your functions in manually at 
# the appropriate point

set.seed(03012018)

source('HW2_sol.R')

# Each of the checks below should return TRUE

# Question 1

# Check 1:

data1 = merge.files('age.txt','teeth.txt')

data1[3,] == c(1,18,28)
is.na(data1[10,3])
is.na(data1[11,2])

# Check 2:
 
data2 = order1(data1,1) 

prod(diff(data2[,1]) > 0) == 1
 
data2[10,] == c(10,18,29)
 
 
# Check 3: 

data3 = order2(data2,2,3)

data3[2,] == c(6,17,28)
data3[5,] == c(10,18,29)

# Question 2


a = seq(0,1,len=15)
b = seq(0,1,len=25)

c = two.sort(a,b)

all.equal(c[27],2/3)
prod( diff(c) >= 0) == 1


# Question 3

# Part a)
# Work with 2(x-3)(x+3)

a2 = c(-18,0,2)
r2 = c(2,3,-3)

standard(2,a2) == -10
horner(2,a2) == -10
factored(2,r2) == -10

# Part b)

all.equal(y.standard[13],-1.818989e-12)
all.equal(y.horner[52],6.252776e-12)
all.equal(y.factored[37],-2.066105e-17)

# Part d)

# Note that these might change by one depending
# on how you counted an exact zero. 

length(zero.standard) == 45
all.equal(zero.horner[25],2.01)
length(zero.factored) == 1


# Question 4

# Part a

set.seed(0623.1912)
datalist = lifedata(20)

all.equal(datalist[[1]][13], 1.1890633, tol=1e-6)

# Part b

# Note that since we have run using a single
# seed from the start of your file, not all
# of these answers will be within tolerance if, 
# for example, you used random numbers in a different
# order.  
#
# Be sensible in assessing whether you have satisfied
# the following. 

all.equal(meandiff, 0.5, tol=1e-2)
all.equal(sd.meandiff, 0.183, tol=1e-2)
all.equal(npts.diff, 33856, tol = 1e1)
all.equal(sd.trunc, 0.011836, tol=1e-3)
all.equal(sd.remove, 0.017414, tol=1e-3)

# Question 5

# We'll try the glider gun configuration

set.seed(1210.1815)

edat = rexp(20)
all.equal(Gini(edat), 0.4432786,tol=1e-6)

# Again, these might be off depending on how 
# you used random numbers

all.equal(gmean, 0.4760605, 1e-6)
all.equal(gvar, 0.00377, 1e-3)

gini.alpha == 0.021

# But this should be exact if you run from the
# last set.seed to here. 

all.equal(gini.power(2, 0.5, 0.06, 1000), 0.006)