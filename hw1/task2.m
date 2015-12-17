function task2
global h
iterMax = 100;
N = 120;
h = 1./N ;
x = 0:h:1;
e = zeros(N+1,1);
errnew = 0;
nn = 1;
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
   
    %     if(err<tol)
    %         break;
    %     end
    e = abs(u-ux); 
%     if(mod(iter,4)==0)
        %        subplot(2,2,nn)
%         figure(1)
%         hold on
%         plot (x, u,'.-b')
%         title(['iter=',num2str(iter)])
%         axis([0 1 0 .5])
%         figure(2)
%         subplot(2,2,nn)
%         semilogy(e)
%         title(['iter=',num2str(iter)])
%         set(gca,'ytick',[1e-15 1e-10 1e-5 1e0])
%         xlabel('x')
%         ylabel('error')
%         axis([0 N 1e-18 1e0])
%         grid on
%         nn = nn + 1;
%         if(nn>4) 
%             break; 
%         end
%    end
err = sum(u-ux)/N;
errnew/err
errnew = err;
end



%         plot (x, u,'.-b')
%         title(['iter=',num2str(iter)])
%         axis([0 1 0 .5])
% iter-1
% 

% figure(2)
% semilogy(e)
% title(['iter=',num2str(iter)])
% set(gca,'ytick',[1e-15 1e-10 1e-5 1e0])
% xlabel('x')
% ylabel('error')
% axis([0 N 1e-18 1e0])
% grid on
end




