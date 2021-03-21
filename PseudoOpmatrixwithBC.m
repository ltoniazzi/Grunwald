function M = PseudoOpmatrixwithBC(BC,fun,n,dx,order)
% This function computes the entries of the (n)x(n) Grunwald
% matrix for the Levy derivative operators with symbol z^2\psi(z) on L_1 [-1,1]
% \psi is the Laplace exponent of the underlying onesided process
% Follows https://arxiv.org/abs/2012.10864
% Example: if \psi(z)=z^\alpha +z^2, 1<\alpha<2, the nonlocal operator is a right Caputo 
% or R-L alpha fractional derivative + a second derivative.
% n denotes the size of the matrix,
% dx
% BC is one of the following cases
% (D^\psi_c,DD) BC=1 (D^\psi_c,DN) BC=2
% (D^\psi_c,ND) BC=3 (D^\psi_c,NN) BC=4
% (D^\psi,ND) BC=5 (D^\psi,NN) BC=6

%% build lower triangular n by n Grunwald matrix
w=LubichWeights(@(s) fun(s),n+1,dx,order);
M=toeplitz(w(2:end),[w(2),w(1),zeros(1,n-2)]);

%% sum(A,2) adds the row entries while sum(A) adds the column entries
% M is the Transition matrix for L1 with default BC=1
%% Boundary weights
switch BC
    case 1
        M(1,:)=0;
        M(end,:)=0;
    case 2 %DalphacDN
        M(end,:)=-sum(M(1:end-1,:)); % boundary weights b^r_i
        M(end,1)=M(end,1)-M(1,2); % overwrite b_n
        M(1,:)=0;
    case 3 %DalphacND
        M(:,1)= -sum((M(:,2:end)),2); % boundary weights b^l_i
        M(end,:)=0;% overwrite b_n      
    case 4 %DalphacNN
        M(:,1)= -sum((M(:,2:end)),2); % boundary weights b^l_i
        M(end,:)=-sum(M(1:end-1,:)); % boundary weights b^r_i
    case 5 %DalphaND
        M(1,1)=M(1,1)+M(1,2); % change only b_l_1
        M(end,:)=0;
    case 6 %DalphaNN
        M(1,1)=M(1,1)+M(1,2); % change only b_l_1
        M(end,:)=-sum(M(1:end-1,:)); % boundary weights b^r_i and b_n
end
