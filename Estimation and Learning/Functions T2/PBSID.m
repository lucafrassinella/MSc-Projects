function [Ai, Bi, Ci, Di, Ki, Sv, X_est] = PBSID(u, y, p, f, n)

% Ensure I/O are row vectors:
u = u(:)';  
y = y(:)';

% Extract parameters:
N = size(y,2);
l = size(y,1);
r = size(u,1); 
m = r+l;

% Construct data vector of inputs-outputs at each time instants:
z = [u; y];
Z_p = zeros(p*m, N-p);
for i = 1:p
    Z_p((i-1)*m+1 : i*m, : ) = z(:, i:N+i-p-1);
end

%%%%% Step #1: VARX estimation (Y = [CKp D] * Zp)

% Assemble data-matrices for regression:
Y = y(:,p+1:N);
U = u(:,p+1:N);
Z_p = [Z_p; U];

% Regression:
CKp_D = Y * pinv(Z_p);

%%%%% Step #2: Construct Lambda_K from the estimated CKp (Markov parameters)
LambdaK = zeros(f*l, p*m);
for i = 1:f
    LambdaK((i-1)*l+1 : i*l, p*m-(p-i+1)*m+1 : p*m) = CKp_D(:, 1:(p-i+1)*m);
end

%%%% Step #3: Perform SVD and estimate future state sequence

[~,S,V] = svd(LambdaK*Z_p(1:p*m, :),'econ');
% Singular Values:
Sv = diag(S)';
X_est = sqrt(S(1:n, 1:n))*V(:,1:n)';

%%%% Step #4: Estimate (A, B, K, C, D) from regression

% remove the window sizes from input and output vector
u = u(:,p+1:p+size(X_est,2));
y = y(:,p+1:p+size(X_est,2));

% CD regression: 
z_CD = vertcat(X_est(:,1:end-1), u(:,1:end-1));
CD = y(:,1:end-1) * pinv(z_CD);

% Estimated innovation error sequence:
E = y - CD*[X_est; u];
    
% ABK regression:
z_ABK = vertcat(X_est(:,1:end-1), u(:,1:end-1), E(:,1:end-1));
ABK = X_est(:,2:end) * pinv(z_ABK);

% Extract outputs:
Ai = ABK(:,1:n);
Bi = ABK(:,n+1:n+r);
Ci = CD(:,1:n);
Di = CD(:,n+1:n+r);
Ki = ABK(:,n+r+1:n+r+l);

end