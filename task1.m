function task1
global h
N = 100;
dt = 0.01;
Tend =1;
tol = 1e-12;
h = 1./N;
assert(dt/h<=1)

x = 0:h:1;

% exact solution 
ux = min(abs(x),abs(x-1))';

u = zeros(N+1,1); 
e = [];

%boundary condition
u(1) = 0;
u(N+1) = 0;
w = u;


for t=dt:dt:Tend
    for j=2:N
        w(j) = u(j) - dt*(HG(u(j-1),u(j),u(j+1)) - 1);
    end
    maxerr = h*sum(abs(u-w));
    if(maxerr<tol)
       fprintf('iter is %d\n',t/dt);
       break;
    end
    u = w;
    err = norm(u-ux);
    
    e = [e err];
    
    figure(1)    
    plot (x, u,'.-b')
    axis([0 1 0 0.5])
    pause (0.05)
%     if(err<tol)
%         break;
%     end
end
figure(2)
semilogy(e)
axis([0 t/dt 1e-16 1e0])

end

function v=HG(um1,u_0,up1)
global h
max1 = max(u_0-um1, 0)^2;
min1 = min(up1-u_0,0)^2;
v = 1./h*sqrt(max(max1,min1));
end

