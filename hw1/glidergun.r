https://powcoder.com
代写代考加微信 powcoder
Assignment Project Exam Help
Add WeChat powcoder
https://powcoder.com
代写代考加微信 powcoder
Assignment Project Exam Help
Add WeChat powcoder
glidergun <- function(n) {
    # initial n*n lattice for a glider gun
    # assumes n >= 40
    A <- matrix(0, n, n)
    A[31,13] <- 1
    A[31,14] <- 1
    A[32,12] <- 1
    A[32,16] <- 1
    A[33,11] <- 1
    A[33,17] <- 1
    A[33,25] <- 1
    A[34,1] <- 1
    A[34,2] <- 1
    A[34,11] <- 1
    A[34,15] <- 1
    A[34,17] <- 1
    A[34,18] <- 1
    A[34,23] <- 1
    A[34,25] <- 1
    A[35,1] <- 1
    A[35,2] <- 1
    A[35,11] <- 1
    A[35,17] <- 1
    A[35,21] <- 1
    A[35,22] <- 1
    A[36,12] <- 1
    A[36,16] <- 1
    A[36,21] <- 1
    A[36,22] <- 1
    A[36,35] <- 1
    A[36,36] <- 1
    A[37,13] <- 1
    A[37,14] <- 1
    A[37,21] <- 1
    A[37,22] <- 1
    A[37,35] <- 1
    A[37,36] <- 1
    A[38,23] <- 1
    A[38,25] <- 1
    A[39,25] <- 1
    return(A)
}
