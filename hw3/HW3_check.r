https://powcoder.com
代写代考加微信 powcoder
Assignment Project Exam Help
Add WeChat powcoder
https://powcoder.com
代写代考加微信 powcoder
Assignment Project Exam Help
Add WeChat powcoder
##############################
##### Homework 3 CheckFile####
##############################

# First we'll call your solutions
# if these are giving you an error, you
# can add your functions in manually at 
# the appropriate point

source('HW3.R')

# Question 1

# Part a)

A = distMat(patchdata[,2])
A[5,7] == 3358

# Part b)

B = dist.cor(patchdata[1:4,1],patchdata[1:4,2])
all.equal(B,0.5917449,1e-4)


# Part c)

set.seed(0217.1890)

C = perm.dist.cor(patchdata[,1],patchdata[,2],10) 

all.equal(C$null.dist[3],0.2174274,1e-4) 
all.equal(C$obs.correlation,0.8907164,1e-4)
  
# Part d) 

set.seed(1810.1919)

D = Q1dfunction(10,10,20)
D == 0.05

# Part e) 

set.seed(0204.1947)

E = Q1efunction(10,10,20)
all.equal( E, c(0.5,0.25) )

# Question 2


# Part a)

set.seed(0105.1949)

A = fdr.data(20,4,4,8)
all.equal(A$p.values[10],0.7989622,tol=1e-4)
all.equal(A$q.values[3],8.673861e-06,tol=1e-4) 

# Part b)

set.seed(0430.2013)
B = fdr.sim(20,4,4,10,10,0.1,0.1)

mean(B$num.genes.selected) == 4.8

all.equal(mean(B$fdp), 0.1361905, 1e-4)
all.equal(sd(B$fdr), 0.08600809, 1e-4)

set.seed(1652.2712)

C = fdr.sim2(20,4,4,10,10,0.1,0.1)
all.equal(C$mean.fdp, 0.1628571,tol=1e-4)
all.equal(C$sd.fdr, 0.07117459,tol=1e-4)

# Part c)

set.seed(3106.2003)
C = Q2cfunction( B,4 )
C$real.one == 1


# Question 3

# Part a)

set.seed(0901.1915)
A = Q3a.function()
all.equal(A[15], 0.2414373, tol=1e-4)

# Part b)

set.seed(2801.1892)
B = Q3simfunc(10)
B$one.significant == 0.8
B$bonferroni.significant == 0



# Question 4

# Part a)

set.seed(2405.1938)
A = ratio.boot(patchdata,10)

all.equal(A$obs.stat, -0.0713061, 1e-4)
all.equal(A$boot.stats[5], -0.148792892, 1e-5)

# Part b)

B = ratio.bias(A)
all.equal(B, 0.001922276, 1e-4)

# Part c)

C = ratio.confint(A)
all.equal(C, c(-0.24606844,0.04217987), 1e-4)

# Question 5

set.seed(1606.1915)
A = PartialPermutation(peanuts,10)

all.equal(A$obs.coef, 0.05580995, 1e-4)
all.equal(A$crit.value, 0.01576853, 1e-4)