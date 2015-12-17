clc; clear all;

N = 50;
a = -1;
L = 2*pi;
dx = L/(N+1);
x = (0:dx:L-dx)';
dt = 2*pi/500;
Tend = 2*pi;
sigma = a*dt/dx;
% CFL 
assert(abs(sigma)<1)

% IC
m1 = (N+1)/4;
m2 = (N+1)/2;
U1 = cos(x)+1/4*sin(m1*x);
U2 = cos(x)+1/4*sin(m2*x);
% plot(x,U2,'b');
alpha = -sigma;

up = .5*alpha*(alpha-1);
do = .5*alpha*(alpha+1);
dia= 1-alpha^2;

e  = ones(N+1,1);
Q = spdiags([do*e dia*e up*e],-1:1,N+1,N+1);

Q(1,end) = do;
Q(end,1) = up;

for t = 0:dt:Tend
    axis([0 2*pi -1 1])
    U1 = Q*U1;
    U2 = Q*U2;
end
hold on
% plot(x,U2,'r*','markersize',3);
plot(x,U2,'go','markersize',3);
axis([0 2*pi -1 1])