% In this code, we solve a 1d wave equation with a 
% jump condition at the interface.
%     u_tt =     u_xx,    xmin<=x<=0, t>=0,
%     w_tt = c^2 u_xx,    0<=x<=xmax, t>=0,
% with initial conditions
%     u(x,0) = U_0(x),  u_t(x,0) = U_1(x),   xmin<=x<=0,
%     w(x,0) = W_0(x),  w_t(x,0) = W_1(x),   0<=x<=xmax,
% and Dirichelt boundary condition at the beginning and Neumann at the end
%     u(xmin,t) = sin(8*pi*t),  w_x(xmax,t) = 0,  t>=0.
% At the interface we have a jump condition
%     u(0,t)   = w(0,t),
%     u_x(0,t) = c^2 w_x(0,t).
% Thomas Frachon, Davoud Saffar 
% 2015-12-29

function wave1d_embed
clf;

%% set the parameters
c    = 2;                 % second wave speed other than zero.
xmin = -1;                % boundary limit
xmax = 1;                 % boundary limit
N    = 56;                % grid points
h    = (xmax-xmin)/(N-1); % space step
dt   = .01;
assert(c^2*dt^2/h^2<1)
Tend =  10;                % final time
filename = sprintf('wave1d_dirichlet_c%d.gif',c*10);

x = linspace(xmin,xmax,N);
% grid points for u and w
Nu = ceil(-xmin/h);
Nw = N-Nu;
% convex parameter for the interpolation at the interface.
alpha = -x(Nu)/h;


% Initial condition
U0 = @(x) 0*ones(numel(x),1);
W0 = @(x) 0*ones(numel(x),1);
U1 = @(x) 0*ones(numel(x),1);
W1 = @(x) 0*ones(numel(x),1);

% v= [u w]';
v  = [U0(x(1:Nu)); W0(x(Nu+1:end))];
v_t= [U1(x(1:Nu)); W1(x(Nu+1:end))];


%% Build two matrices for the system Bv_tt = Av
% Build A
a10 = 1-alpha*(1-c^2);
a11 = -2+(1+alpha)*(1-c^2);
a12 = c^2;
a21 = 1;
a22 = -2+alpha*(1-c^2);
a23 = 1-alpha*(1-c^2);

e1 = ones(N,1);
A  = spdiags([e1 -2*e1 e1], -1:1, N,N)/h^2;
A(Nu,Nu-1:Nu+1) = [a10 a11 a12]/h^2;
A(Nu+1,Nu:Nu+2) = [a21 a22 a23]/h^2;

% Build B
b11 = 1-alpha*(1-c^2) + (alpha-.5)*(1-alpha);
b12 = -(alpha-.5)*(1-alpha);
b21 = alpha*(alpha-.5);
b22 = 1/c^2*(1-alpha*(1-c^2)-c^2*alpha*(alpha-.5));

e1 = ones(Nu,1);
e2 = 1/c^2*ones(Nw,1);
e = [e1; e2];
B   = spdiags(e, 0, N,N);
B(Nu,Nu:Nu+1)   = [b11 b12];
B(Nu+1,Nu:Nu+1) = [b21 b22];


%% time iterations to solve Bv_tt = Av.
% We can convert the system into two 1st order set of ODEs
% as C_1 Z_t = C_2 Z, where Z = [v v_t]'.
% Discretize this by Explicit Euler to get 
% Z^(n+1) = (I+inv(C_1)*C_2*dt*Z^n.
% or Implicit Euler to get 
% Z^(n+1) = inv(I-inv(C_1)*C_2*dt)*Z^n.


Z = [v;v_t];
C = sparse(2*N,2*N);
C(1:N,N+1:2*N)     = eye(N,N);
C(N+1:2*N,1:N)     = B\A;

T = eye(size(C))-dt*C;
T(N,end) = 0;
v = Z(1:N);
v_t = Z(N+1:end);
plot(x,v)

n = 1;i=1;
% time iteration
for t=dt:dt:Tend
    Z = Apply_BC(Z,Nu,Nw,t,h);
    Z = T\Z;
    v = Z(1:N);
    v_t = Z(N+1:end);
    if(mod(n,5)==0)
        plot(x,v)
        hold on
        plot(zeros(N,1),v(Nu),'r.','markersize',20)
        plot(zeros(N,1),linspace(-1,1,N),'k:')
        hold off
        axis([xmin xmax -1 1])
        %drawnow
        frame = getframe();
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if i == 1; 
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append');
        end
        i = i + 1;
    end
    n = n + 1;
end


end

% homogeneous boundary condition
function [UD, UN] = BC(t)
    UD = sin(pi*t);
    UN = 0;
end

function Z = Apply_BC(Z,Nu,Nw,t,h)
    [UD, UN] = BC(t);
    Z(1)     = UD;
    Z(Nu+Nw) = Z(Nu+Nw-1) + h*UN;
end


