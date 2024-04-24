function [U,S,V,R]=new_PbP_QLP_2(A, tol, q) %% m<n
[P] = QB_K(A,tol); % mnd+md^2
if q ~= 0
    for i = 1:q
        C = A'*P;
        [P,~] = qr(C,0);
        C = A*P;
        [P,~] = qr(C,0);
    end
end
D = P'*A; % mnd
[Q,R] = qr(D',0); % nd^2
[P1,R1] = qr(R',0); % d^3
S = R1;
U = P*P1; % md^2
V = Q;
end

%% QB算法
function [Q] = QB_K(A,eps)
[m,n] = size(A);
% Omega = randn(n,m);
Q = zeros(m,n);
R = zeros(n,n);
% count = 5;
% k = round(n/count); %% 设置每次取n/20列进行QR计算
k = 10;
count = round(n/k);
for j = 1:count+1
    ost = (j-1)*k+1; oen = j*k;%min(j*k,n); 
    oenk = ost-1; %% 这三步仅仅是给一个索引
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

% [m,n] = size(A);
% qw = zeros(m,m);
% Omega = randn(n,n);
% v1 = A*Omega(:,1);
% B(1,1) = norm(v1);
% Q(:,1) =v1/B(1,1);
% for j = 2:n
%     for i = 1:j-1
%         sumQ = zeros(size(Q,1));
%         for z = 1:j-1
%             sumQ = sumQ + Q(:,z)*Q(:,z)';
%         end
%         qw(:,j) = A*Omega(:,j)-sumQ*(A*Omega(:,j));
%     end
%     B(j,j) = norm(qw(:,j));
%     if abs(B(j,j)) >= eps
%         Q(:,j) = qw(:,j)/norm(qw(:,j));
%     else
%         break
%     end
% end
end