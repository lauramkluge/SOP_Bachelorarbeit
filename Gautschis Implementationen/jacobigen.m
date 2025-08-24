function J=jacobigen(N,ab)
%   JACOBIGEN(N,AB) generates the NxN array J of
%   recurrence coefficients for the polynomials orthogonal
%   with respect to the weight function W. The
%   weight function W is specified by the Nx2 input array AB
%   of recurrence coefficients for the polynomials orthogonal
%   with respect to the weight function W.

N0=size(ab,1); if N0<N, error('input array ab too short'), end
J=zeros(N);
for n=1:N, J(n,n)=ab(n,1); end
for n=2:N
  J(n,n-1)=sqrt(ab(n,2));
  J(n-1,n)=J(n,n-1);
end