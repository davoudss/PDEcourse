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

cnt = 1;
nn = 1;
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
    
    e = abs(u-ux); 

    if(mod(cnt,12)==0)
         figure(1)
        plot (x, u,'.-b')
        title(['time=',num2str(t)])
        axis([0 1 0 .5])
        figure(2)
        subplot(2,2,nn)
        semilogy(e)
        title(['iter=',num2str(cnt)])
        set(gca,'ytick',[1e-15 1e-10 1e-5 1e0])
        xlabel('x')
        ylabel('error')
        %axis([0 cnt 1e-18 1e0])
        grid on   
        nn = nn + 1;
    end
    cnt = cnt + 1;
 
%     if(err<tol)
%         break;
%     end
end
%         plot (x, u,'.-b')
%         title(['time=',num2str(t)])


end

function v=HG(um1,u_0,up1)
global h
max1 = max(u_0-um1, 0)^2;
min1 = min(up1-u_0,0)^2;
v = 1./h*sqrt(max(max1,min1));
end


