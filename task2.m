function task2
global h
iterMax = 100;
N = 40;
h = 1./N ;
x = 0:h:1;
e = [];
tol = 1e-12;
h = 1./N;

% exact solution
ux = min(abs(x),abs(x-1))';


% initial condition
%u = 2*ones(N+1,1);
Nl = N-1;
Nl1 = ceil(Nl/2);
Nl2 = floor(Nl/2);

u = zeros(N+1,1);
v1 = [ones(Nl1,1); zeros(N-Nl1,1)];
v2 = [zeros(N-Nl2,1) ;ones(Nl2,1)];
B = diag(v1,-1) + diag(v2,1);

%boundary condition
u(1) = 0;
u(N+1) = 0;
b = [u(1)/h; v1] + [v2;u(N+1)/h];

for iter=1:iterMax    
    unew=B*u+h*b;
    err = norm(unew-ux);
    maxerr = h*sum(unew-u);
    if(maxerr<tol)
       fprintf('iter is %d\n',iter);
       break;
    end    
    u = unew;
    
    e = [e err];
%     if(err<tol)
%         break;
%     end
    
    figure(1)
    plot (x, u,'.-b')
    axis([0 1 -2 1])
    
end

iter-1

figure(2)
semilogy(e)
axis([0 iter 1e-12 1e0])
end




