function task3
global h
clf
iterMax = 100;
N = 40;
h = 1./N ;
x = 0:h:1;
Nl = N-1;
Nl1 = ceil(Nl/2);
Nl2 = floor(Nl/2);

u = zeros(N+1,1); 
v0 = [ones(N+1,1)];
v1 = [ones(Nl1,1); zeros(N-Nl1,1)];
v2 = [zeros(N-Nl2,1) ;ones(Nl2,1)];
A1 = diag(v0) + diag(-v1,-1); % + diag(v2,1)
B1 =  diag(v2,1);

v0 = [ones(Nl1+1,1); -ones(Nl2+1,1)];
v1 = [ones(Nl1,1); zeros(N-Nl1,1)];whos v0
v2 = [zeros(N-Nl2,1) ;ones(Nl2,1)];
A2 = diag(v0) + diag(v2,1); % + diag(v2,1)
B2 =  diag(v1,-1);

%boundary condition
u(1) = 0;
u(N+1) = 0;
b1 = [u(1)/h; v1] + [v2;u(N+1)/h];
b2 = [u(1)/h; v1] + [-v2;u(N+1)/h];


for iter=1:iterMax

  utemp = A1\(B1*u+h*b1);
  unew = A2\(B2*utemp+h*b2);
    if(norm(u-unew)<1e-6)
        break;
    else
        u = unew;
    end

    figure(1)    
    plot (x, u,'.-b')
    axis([0 1 -2 1])
    
end
iter-1
end




