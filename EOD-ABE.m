
function [U,S,V,R]=new_PbP_QLP_2(A, tol, q)
[P] = QB_K(A,tol);
if q ~= 0
    for i = 1:q
        C = A'*P;
        [P,~] = qr(C,0);
        C = A*P;
        [P,~] = qr(C,0);
    end
end
D = P'*A;
[Q,R] = qr(D',0);
[P1,R1] = qr(R',0);
S = R1;
U = P*P1;
V = Q;
end

%%
function [Q] = QB_K(A,eps)
[m,n] = size(A);
Q = zeros(m,n);
R = zeros(n,n);
k = 10;
count = round(n/k);
for j = 1:count+1
    ost = (j-1)*k+1; oen = j*k;%min(j*k,n); 
    oenk = ost-1;
    v = A*randn(n,k); % Omega(:,ost:oen);
    v = v - Q(:,1:oenk) * (Q(:,1:oenk)' * v);
    [Q(:,ost:oen),R(ost:oen,ost:oen)] = qr(v,0);
    i = find(abs(diag(R(ost:oen,ost:oen)))<eps,1);
    if ~isempty(i)
        if i == 1 && j ==1
            Q = Q(:,1:k);
            return
        else
            Q = Q(:,1:oen-k+i-1);
            return
        end
    end
end
end
