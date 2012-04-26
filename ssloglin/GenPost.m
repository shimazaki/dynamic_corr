function [MISE,KL_qp,KL_pq,MISEt,KLt_qp,KLt_pq] = GenPost(model,gen)

% Underlying Model
theta1 = gen.param.theta(2:end,:);
eta1 = gen.param.eta(2:end,:);
p1 = gen.param.p;
phi1 = -gen.param.theta(1,:);
psi1 = sum(p1.*log(p1));

% Fitted Model
d = model.binary.d;
N = model.binary.N;
K = model.binary.K;

theta2 = zeros(2^N-1,K); 
eta2 = zeros(2^N-1,K); 
p2 = zeros(2^N,K); 

theta2(1:d,:) = model.param.smoother.theta;
eta2(1:d,:) = model.param.smoother.eta;
p2 = model.param.smoother.p;
phi2 = model.param.smoother.phi;
psi2 = sum(p2.*log(p2));

% Distances, MISE and KL
D = model.struct.D;
MISEt = 0; 
M = length(theta2(:,1));
for i = 1: M
    MISEt = MISEt + 1/M* ( theta1(i,:)-theta2(i,:) ).^2 ;
end
MISE = sum(MISEt)*D;

KLt_qp = phi2 + psi1 - sum(theta2.*eta1);
KL_qp = sum(KLt_qp);

KLt_pq = phi1 + psi2 - sum(theta1.*eta2);
KL_pq = sum(KLt_pq);

%[KL_qp sum( sum(p1.*log(p1./p2) ) )]
%[KL_pq sum( sum(p2.*log(p2./p1) ) )]
%[MISE,KL_qp,KL_pq]

