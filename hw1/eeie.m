clc; clear all;close all

Lambda = [10 90 120 140 190 200];

for i=1:numel(Lambda)
    lambda= Lambda(i);
    h = .01;
    h*lambda
    n = 1:100;
    
    EX = exp(lambda*n*h);
    EE = (1+h*lambda).^(n);
    IE = (1-h*lambda).^(-(n));
    
    EE_err = abs((EX-EE)./EX);
    IE_err = abs((EX-IE)./EX);
    
    subplot(3,2,i)
    semilogy(n,EE_err, '.-b',n,IE_err,'o-r','markersize',3)
    title(['h\lambda=',num2str(h*lambda)])
    legend('EE','IE')
    xlabel('iter')
    ylabel('Relative error')
    % hold on
    % legend('Explicit','Imiplicit')
    % semilogy(n,EE_err-IE_err, '.-b')
end

return
figure(2)
plot(EX,'g')
hold on
plot(EE,'.-b','markersize',3)
plot(IE,'o-r','markersize',3)
xlabel('x')
ylabel('y')
legend('Exact','EE','IE')
