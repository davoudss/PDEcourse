clc; clear all;clf

lambda = 3;
h = .1;

n = 1:10;

EX = exp(lambda*n*h);
EE = (1+h*lambda).^(n+1);
IE = (1-h*lambda).^(-(n+1));

EE_err = (EX-EE)./EX;
IE_err = (EX-IE)./EX;

%semilogy(n,EE_err, '.-b',n,IE_err,'o-r')
%legend('Explicit','Imiplicit')
%semilogy(n,EE_err-IE_err, '.-b')

plot(EX)
hold on
plot(EE,'.-')
plot(IE,'o-g')
legend('Exact','EE','IE')
