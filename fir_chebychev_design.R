# This program designs an FIR filter, given a desired frequency response Hdes(ω). 
# The design is judged by the maximum absolute error (Chebychev norm). 

# A derivative work by Judson Wilson, 5/27/2014 in CVXPY. Adapted from the MATLAB's CVX example of the same name, by Almir Mutapcic, 2/2/2006.
# and now ported to R 

# See http://www.cvxpy.org/en/1.0/examples/applications/fir_chebychev_design.html

n <- 20
j <- sqrt(as.complex(-1))
# rule-of-thumb frequency discretization (Cheney's Approx. Theory book)
m <- 15 * n
w <- as.matrix(seq(0, pi, length.out = m))
D <- 8.25
Hdes = exp(-j*D*w) # desired frequency response
#var = 0.05;
#Hdes = 1/(sqrt(2*pi*var))*exp(-(w-pi/2)^2/(2*var));
#Hdes = Hdes*exp(-j*n/2*w);

A <- t(exp( -j*kronecker(t(w), 1:(n)) ))

# optimal Chebyshev filter formulation
h <-   Variable(n,1)

Hdes_real <- Re(Hdes)
Hdes_imag <- Im(Hdes)

obj <- Minimize( max(
  square(Re(A) %*% h - Hdes_real) + square(Im(A) %*% h - Hdes_imag)
  )  )
soln <- solve(Problem(obj))

h <- soln$getValue(h)

H <- t(exp(-j*kronecker(t(w),1:(n)))) %*% h
plot(w, 20*log10(abs(H)), col = "red", lwd=2,type = "l", ylim = c(-30, 10))
lines(w, 20*log10(abs(Hdes)), col = "blue", lwd=2, type = "l", lty=2, ylim = c(-30, 10))
grid(5,5)
plot(w, Arg(H), type = "l")
