% SOBZEROS Zeros of Sobolev orthogonal polynomials.
%
%    Given the NxN matrix B=B_N of the Sobolev recurrence coefficients, 
%    z=SOBZEROS(n,N,B) generates the (column) vector z of the zeros of 
%    the nth-degree Sobolev orthogonal polynomial, where 1<=n<=N. 
%    The zeros, if all real, are ordered increasingly, otherwise 
%    with increasing moduli.
%
function [v, z] = sobzeros_mine2(n, N, B, scaling)
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
