https://powcoder.com
代写代考加微信 powcoder
Assignment Project Exam Help
Add WeChat powcoder
https://powcoder.com
代写代考加微信 powcoder
Assignment Project Exam Help
Add WeChat powcoder
##############################yu
##### Homework 3 CheckFile####
##############################

# First we'll call your solutions
# if these are giving you an error, you
# can add your functions in manually at 
# the appropriate point

source('HW4_sol.R')

# Question 1

# Part a)

all.equal(opt.x.gs[3],0.1524931,tol=1e-7)
niter.gs[5] == 33

# Part b)

all.equal(opt.x.nr[19],12.3428886,tol=1e-7)
niter.nr[7] == 8

# Part c)

niter.gs.log[5] == 33
all.equal(opt.x.nr.log[15],1.967027,tol=1e-7)

# Part d)

all.equal(opt.x.gs2[1],1.715918,tol=1e-7)
niter.nr2[6] == 9

all.equal(opt.x.nr3[11],1.113138e-01,tol=1e-7)
niter.gs3[9] == 34


# Question 2

# Part a

try1 = muon.golden.section(muon)
all.equal(try1$value, -62.95529, tol=1e-7)

# Part b

set.seed(041318)
all.equal( boot.muon.ci(muon)$confint[1], 0.3083642, tol=1e-6)

# Part c

all.equal( likelihood.ratio.ci(muon,try1$opt.alpha)[2],0.9334721,tol=1e-7)


# Question 3

# Note, could be very sensitive to order of random numbers

# Part a
#
# Note: this was done using the fewest random numbers that 
# Giles could think of. 

set.seed(031918)
A = rtruncnorm1(100,1)
all.equal(A$X[15],1.839489,tol=1e-6)  
A$N == 317


# Part c

set.seed(032018)
B = rtruncnorm2(100,1)
B$N == 140
all.equal(B$X[78], 1.391446,tol=1e-6)


# Part d

# Note 10% tolerance, might still be low

all.equal(effort[3,],c(376327,1088,370398.347),tol=0.1)



# Question 4


# Part a

PoissonReg1(0,0,Flu)$num.iter == 27

all.equal(PoissonReg1(0,0,Flu)$b0,1.990234,tol=1e-6)


# Part b

all.equal(PoissonReg2(0,0,Flu)$likelihood,-194.5541)

all.equal(PoissonReg2(0,0,Flu)$b1, -0.01746317, tol=1e-4)

# Part b

# Very loose tolerance should be ok
all.equal(Flu.confint.b1,c(-0.02077101,-0.01415533), tol=1e-1)