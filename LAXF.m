clc; clear all; close all

N = 100;
a = -1;
L = 2*pi;
dx = L/(N+1);
x = (0:dx:L-dx)';
dt = .01;
Tend = 2;
sigma = a*dt/dx;
% CFL 
assert(abs(sigma)<1)

% IC
m1 = (N+1)/4;
m2 = (N+1)/2;
U1 = cos(x)+1/4*sin(m1*x);
U2 = cos(x)+1/4*sin(m2*x);

up = .5*(1-sigma);
do = .5*(1+sigma);

e  = ones(N+1,1);
Q = spdiags([do*e up*e],[-1 1],N+1,N+1);

Q(1,end) = do;
Q(end,1) = up;

for t = 0:dt:Tend
    plot(x,U1,x,U2,'r');
    axis([0 2*pi -1 1])
    pause(.05)
    U1 = Q*U1;
    U2 = Q*U2;
end

plot(x,U1,x,U2,'r');
axis([0 2*pi -1 1])