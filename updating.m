function H = updating(Z,w,proc)
% Updating procedure to generate polynomials orthogonal with respect to a
% discretized diagonal Sobolev inner product given matrix Z and weigths w
%INPUT
%   Z = matrix in Jordan canonical form
%   w = weight vector
%   proc = 'PR' for plane rotations or 'HH' for Householder reflectors
%OUTPUT
%   H = Hessenberg matrix containing recurrence coefficients for the
%   sequence of Sobolev orthonormal polynomials

m = length(w);
%w = w/norm(w); % war vorher auskommentiert
%% First block
k1 = find(w>0,1); % size of first block
J1 = Z(1:k1,1:k1);
w1 = w(1:k1);

% solution to initial HIEP
H = J1';

%% Iterate and merge all other blocks
% Remaining weights and nodes
wrem = w(k1+1:m);
Jrem = Z(k1+1:m,k1+1:m);
while ~isempty(wrem)
    kl = find(wrem>0,1); % size of first block
    Jl = Jrem(1:kl,1:kl);
    wl = wrem(1:kl);
    wdot = w(1:(m-length(wrem)));

    Jrem = Jrem(kl+1:end,kl+1:end);
    wrem = wrem(kl+1:end);

    % Merge
    if proc=='PR'
        H = merge_PR(H,Jl',wdot,wl);
    elseif proc=='HH'
        H = merge_HH(H,Jl',wdot,wl);
    end
    
end
end

function [Htil]= merge_PR(H1,H2,w1,w2)
% Merges two solutions to a Hessenberg IEP to form solution to a larger Hess IEP
% INPUT:
%       H1,H2 = Hessenberg matrix of solution of some Hess IEP
%       w1,w2 = (original) weigths of the Hess IEPs
% OUTPUT:
%       Htil = solution to Hess IEP obtained by combination of two given IEPs

n = length(w1);
m = length(w2);
I = eye(n+m);

% Merge the given matrices and weigts
Htil = [H1, zeros(n,m);...
        zeros(m,n),H2];
wtil = [w1;w2]; 
   
% Enforce orthogonality condition with respect to new weigth vector
P= I; P([1,n+1],[1,n+1]) = givens(norm(w1),norm(w2))';
Htil = P'*Htil*P;

%% Restore Hessenberg structure column by column
for i=2:n
  for j=min(m,i):-1:2
    P= I; P([n+j-1,n+j],[n+j-1,n+j]) = givens(Htil(n+j-1,i-1),Htil(n+j,i-1))';
    
    Htil = P'*Htil*P;  
    Htil(n+j,i-1) = 0; % set explicitly to zero
  end
  P= I; P([i,n+1],[i,n+1]) = givens(Htil(i,i-1),Htil(n+1,i-1))';
  Htil = P'*Htil*P;
  Htil(n+1,i-1) = 0; % set explicitly to zero
end
i = 1;
for j = m:-1:i+1
    P= I; P([n+j-1,n+j],[n+j-1,n+j]) = givens(Htil(n+j-1,n+i-1),Htil(n+j,n+i-1))';
    Htil = P'*Htil*P;
    Htil(n+j,n+i-1) = 0; % set explicitly to zero
end
for i=2:m-1
  for j = m:-1:i+1
    P= I; P([n+j-1,n+j],[n+j-1,n+j]) = givens(Htil(n+j-1,n+i-1),Htil(n+j,n+i-1))';
    Htil = P'*Htil*P;
    Htil(n+j,n+i-1) = 0; % set explicitly to zero
  end

end
end

function [Htil]= merge_HH(H1,H2,w1,w2)
% Merges two solutions to a Hessenberg IEP to form solution to a larger Hess IEP
% To restore the structure Householder reflectors are used
% INPUT:
%       H1,H2 = Hessenberg matrix of solution of some Hess IEP
%       w1,w2 = (original) weigths of the Hess IEPs
% OUTPUT:
%       Qtil = solution to Hess IEP obtained by combination of two given IEPs

n = length(w1);
m = length(w2);
I = eye(n+m);

% Merge the given matrices and weigts
Htil = [H1, zeros(n,m);...
        zeros(m,n),H2];
wtil = [w1;w2]; 
   
% Enforce orthogonality condition with respect to new weigth vector
P= I; P([1,n+1],[1,n+1]) = givens(norm(w1),norm(w2))';
Htil = P'*Htil*P;

%% Restore Hessenberg structure column by column
for i=1:n-1
  irow = min(i+1,m); 
  x = [Htil(i+1,i);Htil([n+1:n+irow],i)];
  omega = angle(x(1));
  e = exp(1i*omega);
  v = (x-e*norm(x)*[1;zeros(irow,1)])/(norm(x-e*norm(x)*[1;zeros(irow,1)]));
  Pdot = eye(irow+1)-2*v*v';
  P = eye(n+m);
  P([i+1,n+1:n+irow],i+1) = Pdot(:,1);
  P([i+1,n+1:n+irow],[n+1:n+irow]) = Pdot(:,2:irow+1);
  Htil = P*Htil*P';
  Htil([n+1:n+irow],i) = 0; % set explicitly to zero
end
ifinal = n+m;
i = n;
while i<n+m-1
  if n+i+1<=n+m
    ifinal = n+i+1;
  end 
  x = Htil([i+1:ifinal],i);
  omega = angle(x(1));
  e = exp(1i*omega);
  v = (x-e*norm(x)*[1;zeros(ifinal-i-1,1)])/(norm(x-e*norm(x)*[1;zeros(ifinal-i-1,1)]));
  Pdot = eye(ifinal-i)-2*v*v';
  P = eye(n+m);
  P([i+1:ifinal],[i+1:ifinal]) = Pdot;
  Htil = P*Htil*P';
  Htil([i+2:ifinal],i) = 0; % set explicitly to zero
  
  i = i + 1;
end
end
  
function G=givens(x,y)
% Computes Givens rotation P such that P [x;y] = [*;0], with * a nonzero element
if x==0
    c=0;
    s=1;
    G=[c,s;-conj(s),c]; % von mir hinzugefügt
    return;
end

aa=abs(x);
scale=aa+abs(y);
nrm=scale*sqrt(abs(x/scale)^2+abs(y/scale)^2);
alpha=x/aa;
c=aa/nrm;
s=alpha*conj(y)/nrm;

G=[c,s;-conj(s),c];
end
