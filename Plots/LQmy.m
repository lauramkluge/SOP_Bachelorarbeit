% LEASTSQAURESAPPROX Computes a least squares approximation to a given
% function using the Legendre measure and the Althammer measure.
% This code generates Figure 5.
close all
clearvars

m=201; g=1/100;
ab=r_jacobi(m); 
zw=gauss(m,ab);

%% Arnoldi
% The Arnoldi iteration is used throughout this experiment
% Starting vector
wsob = zeros(2*m,1); wsob(2:2:2*m) = sqrt(zw(:,2));
% Jordan matrix
Zsob = zeros(2*m);
Zsob(1:2:end,1:2:end) = diag(zw(:,1));
Zsob(2:2:end,2:2:end) = diag(zw(:,1));
for k=1:2*m
    if mod(k,2)==1
       Zsob(k,k+1) = sqrt(g);
    end
end

% Compute recurrence matrix and orthonormal basis using Arnoldi
[Vsob,Hsob] = Arnoldi(Zsob,wsob,m); % For the Althammer measure
HsobUp = updating(Zsob, wsob, 'PR');
size(Vsob)
[n,m] = size(Hsob);
VsobUp = sylvester(Zsob, -HsobUp(1:202, 1:202), zeros(402, 202));

%[V,H] = Arnoldi(diag(zw(:,1)),sqrt(zw(:,2)),m); % For the Legendre measure


%% Solve least squares problem
omega = 100;
f = @(x) exp(-omega*(x-1/5).^2); % Function of interest (FOI)
der = @(x) -2*omega*(x-1/5).*exp(-omega*(x-1/5).^2); % Derivative of FOI

ff = f(zw(:,1));
derder = der(zw(:,1));
rhssob = zeros(2*length(zw),1); rhssob(1:2:end) = sqrt(g)*sqrt(zw(:,2)).*derder; rhssob(2:2:end) = sqrt(zw(:,2)).*ff;
rhs = sqrt(zw(:,2)).*ff;

%% Plotting error of the least squares approximants to the function and its derivative
p0 = 1/norm(wsob);

xx = linspace(-1,1,1000);

err_fctsob = [];
err_dersob = [];
err_fct = [];
err_der = [];
nn = 1:10:m;
for ind=1:length(nn)
    i = nn(ind);
    Qsob = Vsob(:,1:i);
    Q = VsobUp(:,1:i);
    
    csob = Qsob'*rhssob;
    c = Q'*rhssob;

    evalsob = polyvalJ(Hsob(1:i,1:i-1),xx,xx,p0,csob);
    eval = polyvalJ(HsobUp(1:i,1:i-1),xx,xx,p0,c);

    err_fctsob = [err_fctsob, norm(f(xx(:))-evalsob(2:2:end),"inf")/norm(f(xx),"inf")];
    err_fct = [err_fct, norm(f(xx(:))-eval(2:2:end),"inf")/norm(f(xx),"inf")];
    
    err_dersob = [err_dersob, norm(der(xx(:))-evalsob(1:2:end),"inf")/norm(der(xx),"inf")];
    err_der = [err_der, norm(der(xx(:))-eval(1:2:end),"inf")/norm(der(xx),"inf")];

end

fig = figure
subplot(1,2,1)
semilogy(nn,err_fctsob,'b*')
hold on
semilogy(nn,err_fct,'ro')
xlabel('n')
ylabel('err in infnorm')
legend hide
subplot(1,2,2)
semilogy(nn,err_dersob,'b*')
hold on
semilogy(nn,err_der,'ro')
xlabel('n')
ylabel('err in infnorm')
legend hide