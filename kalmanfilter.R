
# This is a robust Kalman Filter Example using CVXR. This is the initial dynamic simulation. Another script will try to use CVXR to solve
# the standard Kalman filtering problem to recover the system's states (x_true)
# The other script is kalmanfilter2.R
#
# This program is ported from CVXPY by Steve Diamond
# See http://nbviewer.jupyter.org/github/cvxgrp/cvxpy/blob/master/examples/notebooks/WWW/robust_kalman.ipynb


n = 1000 # number of timesteps
T1 = 50 # time will vary from 0 to T with step delt
ts = pracma::linspace(0, T1, n)
gamma = .05 # damping, 0 is no damping
delt = ts[2]-ts[1]


A = pracma::zeros(4,4)
B = pracma::zeros(4,2)
C = pracma::zeros(2,4)

A[1,1] = 1
A[2,2] = 1
A[1,3] = (1-gamma*delt/2)*delt
A[2,4] = (1-gamma*delt/2)*delt
A[3,3] = 1 - gamma*delt
A[4,4] = 1 - gamma*delt

B[1,1] = delt^2/2
B[2,2] = delt^2/2
B[3,1] = delt
B[4,2] = delt

C[1,1] = 1
C[2,2] = 1

sigma = 20
p = .20
set.seed(6)

x = pracma::zeros(4,n+1)
x[,1] = c(0,0,0,0)
y = pracma::zeros(2,n)

# generate random input and noise vectors
w = pracma::randn(2,n)
v = pracma::randn(2,n)

# add outliers to v
set.seed(0)
inds = pracma::rand(n) <= p
inds = inds*1
v[,inds] = sigma*pracma::randn(2,n)[,inds]

# simulate the system forward in time
for (t in 1:n){
  y[,t] = C%*%(x[,t]) + v[,t]
x[,t+1] = A%*%(x[,t]) + B%*%(w[,t])
}
x_true = x
w_true = w
#plot(ts, x_true[1,-1], type = "l")
# plot(ts, x_true[2,-1], type = "l")
# plot(ts, x_true[3,-1], type = "l")
# plot(ts, x_true[4,-1], type = "l")
 plot(ts, w_true[1,], type = "l")
# plot(ts, w_true[2,], type = "l")
