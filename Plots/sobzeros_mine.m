% SOBZEROS Zeros of Sobolev orthogonal polynomials.
%
% Entspricht der Funktion sobzeros.m (siehe Ordner: Gautschis
% Implementationen) mit der Änderung, dass eine Skalierung der Matrix zur
% besseren Konditionierung des EW-Problems genutzt werden kann. Außerdem
% wird zusätzlich noch die Basis aus Eigenvektoren zurückgegeben.
function [v, z] = sobzeros_mine(n, N, B, scaling)
if(n<1|n>N), error('n out of range'), end
H=zeros(n);
for i=1:n
  for j=1:n
    if i==1
      H(i,j)=B(j,j);
    elseif j==i-1
      H(i,j)=1;
    elseif j>=i
      H(i,j)=B(j-i+1,j);
    end
  end
end
if nargin<4
    scaling = ones(n,1);
end
[v,z]=eig(diag(scaling) * H * diag(1./scaling), "vector");