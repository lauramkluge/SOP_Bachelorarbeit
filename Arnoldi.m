function [V,H] = Arnoldi(A,v,n)
%ARNOLDI Arnoldi iteration
%   Generates an orthonormal basis V for the Krylov subspace 
%   span{v,Av,A^2 v,...}
%   It also generates the Hessenberg recurrence matrix H = V' A V
%INPUT:
%   A = mxm matrix
%   v = column vector of size m
%OUTPUT
%   V = orthonormal basis for span{v,Av,A^2 v,...,A^(m-1)v}
%   H = Hessenberg matrix satisfying A*V=V*H
m = length(A);

v = v(:); % make column vector
v = v/norm(v);

V = [v];
if n<m
    H = zeros(n+1,n);
else
    H = zeros(m);
end
for i=1:n
    vhat = A*V(:,i);
    for j = i:-1:1
        h = dot(V(:,j),vhat);
        vhat = vhat-h*V(:,j);
        H(j,i) = h;
    end
    if i<m
        H(i+1,i) = norm(vhat);
        V = [V,vhat/norm(vhat)];
    end
end

end
