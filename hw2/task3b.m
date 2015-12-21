function task3b
clc; clear all;

% input data
k = 3;
nn = [1 4 7];
NN = [40 80 120];
n = nn(k);
N = NN(k);

scheme = '3-point';

% plot vectors
x = linspace(-.5,.5,N);
h = x(2)-x(1);
y = x;
[x,y] = meshgrid(x,y);

subplot(3,3,n)
uex = func(x,y);

mesh(x,y,uex)
colormap jet
view(2)

% Matrix construction
if(strcmp(scheme,'3-point'))
    A = Construct_Matrix_scheme2(N);
    % Apply BC
    b = Apply_BC_scheme2(N,x,y);
else
    A = Construct_Matrix_scheme1(N);
    % Apply BC
    b = Apply_BC_scheme1(N,x,y);
end
if(N<10)
    full(A)
    full(b)
end
b = reshape(b,N^2,1);
u = A\b;
u = reshape(u,N,N);
subplot(3,3,n+1)
mesh(x,y,u)
colormap jet
view(2)

subplot(3,3,n+2)
mesh(x,y,max(log(abs(u-uex)),-16))
colormap jet
colorbar
view(2)
end

function A = Construct_Matrix_scheme1(N)
A = zeros(N^2,N^2);
A= sparse(A);
e = ones(N,1);
for k=1:N:N^2
    if(k==1 || k==N^2-N+1)
        A(k:k+N-1,k:k+N-1) = diag(e);
    else
        A(k,k) = 1; A(k+N-1,k+N-1)=1;
        for j=1:N-2
            A(k+j,k+j) = -6;     %diagonal
            A(k+j,k+j+1) = 2;    %right
            A(k+j,k+j-1) = 2;    %left
            A(k+j,k+j+N) = 2;    %up
            A(k+j,k+j-N) = 2;    %down
            A(k+j,k+j+N-1) = -1; %up-left
            A(k+j,k+j-N+1) = -1; %down-right
        end
    end
    
end
end

function A = Construct_Matrix_scheme2(N)
A = zeros(N^2,N^2);
A= sparse(A);
e = ones(N,1);
for k=1:N:N^2
    if(k==1 || k==N^2-N+1)
        A(k:k+N-1,k:k+N-1) = diag(e);
    else
        A(k,k) = 1; A(k+N-1,k+N-1)=1;
        for j=1:N-2
            A(k+j,k+j) = -2;    %diagonal
            A(k+j,k+j+N+1) = 1; %up-right
            A(k+j,k+j-N-1) = 1; %down-left
        end
    end
    
end
end


function b=Apply_BC_scheme1(N,x,y)
b = zeros(N,N);
v = 1:N:N^2;
b(v)=func(x(v),y(1));
v = 1+N-1:N:N^2;
b(v)=func(x(v),y(end));
v = 1:1+N-1;
b(v) = func(x(1),y(v));
v = N*(N-1)+1:N^2;
b(v) = func(x(end),y(v));
end

function b=Apply_BC_scheme2(N,x,y)
b = zeros(N,N);
v = 1:N:N^2;
b(v)=func(x(v),y(1));
v = 1+N-1:N:N^2;
b(v)=func(x(v),y(end));
v = 2:1+N-2;
b(v) = func(x(1),y(v));
v = N*(N-1)+2:N^2-1;
b(v) = func(x(end),y(v));
end

function f = func(x,y)
a = 0;
f = (x<=y).*sin(6*pi*(x-y)) + a*(x>y).*sin(6*pi*(x-y)) ;
end