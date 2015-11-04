function task2
global h
iterMax = 20;
N = 20;
h = 1./N ;
x = 0:h:1;
% initial condition
%u = 2*ones(N+1,1);
Nl = (N+1)/2; 
u = zeros(N+1,1); 
v1 = [ones(floor(Nl),1); zeros(N-floor(Nl),1)];
v2 = [zeros(N-floor(Nl)+1,1) ;ones(floor(Nl)-1,1)];
B = diag(v1,-1) + diag(v2,1);

%boundary condition
u(1) = -2;
u(N+1) = -2;
b = [u(1)/h; v1] + [v2;u(N+1)/h];
for iter=1:iterMax

        u=B*u+h*b;

    figure(1)    
    plot (x, u,'.-b')
    axis([0 1 -2 1])
    
    pause (0.1)
end


end




