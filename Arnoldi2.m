function [Q,H] = Arnoldi2(A,b,k)

n = length(A);
Q = zeros(n,k+1); % Orthonormal basis
H = zeros(k+1,k); % Upper Hessenberg matrix
alpha = 0.1; % The "twice is enough" threshold

Q(:,1) = b/norm(b);
for j = 1:k

% Get a vector in the next subspace (and its norm)
Q(:,j+1) = A*Q(:,j);
norma = norm(Q(:,j+1));

% Modified Gram-Schmidt (standard Arnoldi)
for l = 1:j
H(l,j) = Q(:,l)'*Q(:,j+1);
Q(:,j+1) = Q(:,j+1)-Q(:,l)*H(l,j);
end
H(j+1,j) = norm(Q(:,j+1));

% The "twice is enough" second pass, if the residual is small
if H(j+1,j) < alpha*norma
for l = 1:j
mu = Q(:,l)'*Q(:,j+1);
Q(:,j+1) = Q(:,j+1)-Q(:,l)*mu;
H(j,l) = H(j,l) + mu;
end
H(j+1,j) = norm(Q(:,j+1));
end

% Normalize final result
Q(:,j+1) = Q(:,j+1)/H(j+1,j);

end
end