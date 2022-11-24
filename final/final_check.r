https://powcoder.com
代写代考加微信 powcoder
Assignment Project Exam Help
Add WeChat powcoder
https://powcoder.com
代写代考加微信 powcoder
Assignment Project Exam Help
Add WeChat powcoder
##############################
##### Final Exam CheckFile####
##############################

# First we'll call your solutions
# if these are giving you an error, you
# can add your functions in manually at 
# the appropriate point

source('final_sol.R')

# Question 1


# Part a

set.seed(11011807)
X = rnorm(100)

eq = function(x){x}
all.equal(MCfunc1(X,sin,eq),c(0.0041542073, -0.0003166141),tol=1e-7)


# Part d

set.seed(11071832)

testfn = function(x){ cos(x) }
X = rnorm(7)

all.equal( as.vector(hfunc(c(0,0.2),X,0.3,testfn)), c(1.013622,1.116312), tol=1e-7)


# Question 2

# Part a
set.seed(12071844)

X = rnorm(10)
all.equal(as.vector(kde(c(-0.2,1),X,0.5)),c(0.3772147,0.2583243),tol=1e-6)

# Part b

set.seed(12071844)

X = rnorm(10)
all.equal(kde.cv(X,0.34),-12.44386,tol=1e-6)


# Part c

# Note use of random numbers here; I chose which observations
# first, and then simulated normals around them. 

set.seed(12071844)

X = rnorm(10)
all.equal(kde.sim(2,X,0.45),c(0.04231365,-0.70325134))

# Part d

set.seed(12071844)

X = rnorm(10)
Xsamp = rnorm(10)

all.equal(Afn(0.2,X,Xsamp,0.52),1.106754,tol=1e-6)

## No examples from maximization or MCMC: new samples
## are drawn each value of theta. 


# Question 3

# Part a
all.equal(logistic.nll(c(-2.5,-1.5,1.1),toenail),142.258,tol=1e-6)

# Part b
set.seed(2106.1864)
all.equal(logistic.nll.mc(c(-2.5,-1.5,1.1),toenail),140.326,tol=1e-6)

# Part d

set.seed(2106.1864)
Z = rnorm(12,sd=0.5)
all.equal(joint.posterior(c(0.1,-0.25),Z,toenail),-69.32814,tol=1e-7)

# Part e

# Usual warning about random numbers inside your functions.
# Although less likely to occur here. 

set.seed(2106.1864)
all.equal(ranef.post(1000,c(-0.19,-0.21,1),Subj1)[1:3],
      c(1.07871834,0.12607701,0.05866123),tol=1e-6)
