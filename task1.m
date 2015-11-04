function task1
global h
N = 20;
dt = 0.1;
Tend =1;
h = 1./N ;
x = 0:h:1;
u = zeros(N+1,1); 

%boundary condition
u(1) = 0;
u(N+1) = 0;
w = u;

for t=dt:dt:Tend
    for j=2:N
        w(j) = u(j) - dt*(HG(u(j-1),u(j),u(j+1)) - 1);
    end
    u = w;

    figure(1)    
    plot (x, u,'.-b')
    axis([0 1 0 0.5])
    pause (0.05)
end


end

function v=HG(um1,u_0,up1)
global h
max1 = max(u_0-um1, 0)^2;
min1 = min(up1-u_0,0)^2;
v = 1./h*sqrt(max(max1,min1));
end


