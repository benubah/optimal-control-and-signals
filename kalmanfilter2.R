# http://nbviewer.jupyter.org/github/cvxgrp/cvxpy/blob/master/examples/notebooks/WWW/robust_kalman.ipynb

# Run kalmanfilter.R first to obtain the values for y

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

B[1,1] = (delt^2)/2
B[2,2] = (delt^2)/2
B[3,1] = delt
B[4,2] = delt

C[1,1] = 1
C[2,2] = 1

#y = pracma::zeros(2,n)

x = Variable(4,n+1)
w = Variable(2,n)
v = Variable(2,n)

tau = .08

obj = sum_squares(w) + tau*sum_squares(v)
obj = Minimize(obj)
constr = c()
for (t in 1:n){
  constr1 = list(x[,t+1] == A%*%x[,t] + B%*%w[,t],
                y[,t]   == C%*%x[,t] + v[,t])
  constr = c(constr, constr1)
 }
prob = Problem(obj, constr)
sol = solve(prob, verbose=TRUE)
x1 = sol$getValue(x)
w1 = sol$getValue(w)
plot(ts, x1[1,-1], type = "l", col = "blue") # optimal x-position
lines(ts, x_true[1,-1], type = "l", col = "red") # true x-position
plot(ts, x1[2,-1], type = "l", col = "blue") # optimal y-position
lines(ts, x_true[2,-1], type = "l", col = "red") # true y-position
plot(ts, x1[4,-1], type = "l", col = "red") # optimal y-velocity
lines(ts, x_true[4,-1], type = "l", col = "blue") # true y-velocity
plot(ts, w_true[1,], type = "l", col = "blue", ylim = c(-3,3))
lines(ts, w1[1,], type = "l", col = "red")

