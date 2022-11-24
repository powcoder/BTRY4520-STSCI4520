https://powcoder.com
代写代考加微信 powcoder
Assignment Project Exam Help
Add WeChat powcoder
https://powcoder.com
代写代考加微信 powcoder
Assignment Project Exam Help
Add WeChat powcoder
##############################
##### Homework 1 CheckFile####
##############################

# First we'll call your solutions
# if these are giving you an error, you
# can add your functions in manually at 
# the appropriate point

source('HW1.R')

# Each of the checks below should return TRUE

# Question 1

# Check 1:

all.equal(y.values[13], 0.512) 

# Check 2:

all.equal(fn(1/sqrt(2)),0.5)

# Check 3: 
all.equal(fn(1.3),sqrt(1.3))

# Question 2

all.equal(hn.for(0.5,6),1.984375)
all.equal(hn.while(0.5,6),1.984375)
all.equal(hn.vec(0.5,6),1.984375)
all.equal(hn.skip(0.5,6),1.328125)

# Question 3

A = matrix(1:6,3,2)
B = matrix(1:4,2,2)

all.equal(matprod(A,B),A%*%B)


# Question 4

growth = read.csv('treegrowth.csv')

# You can check your solution with the plot in the book.

# But we can check an extracted slope

all.equal(as.numeric(ageslope(15,growth)), 1.128) 


# Question 5

# We'll try the glider gun configuration

source('glidergun.r')

A  = glidergun(50)

# Checks on specific values

neighbours(A,35,35,50) == 2
neighbours(A,1,50,50) == 0
neighbours(A,35,34,50) == 1