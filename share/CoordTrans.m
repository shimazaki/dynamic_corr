function [th, eta, p, rho] = CoordTrans(theta,N,R,Pr,THr,invTHr)
% This function convenverts 
% input
% theta: a log-linear parameters
% N: number of binary variables
% R: the order of interactions
% Pr, THr, invTHr: matrices for coordinate transformation
% output
% th: a log-lienar parameters with the first row being the normalization parameter.
% eta: eta-coordinate, probabilites of synchornies.
% p: p-coordinate, probabilities of spike patterns.
% rho: marginal correlations
M = 0; 
for r = 1: R
    M = M + nchoosek(N,r);
end

if ( nargin > 3 )
    PR = Pr{N}(1:M+1,:);
    THR = THr{N}(:,1:M+1);
else
	for r = 1 : R
		[THr{r},IDX{r}]= CoordGen(r);
		Pr{r} = THr{r}';
		invTHr{r}=inv(THr{r});
	end
	
	[THR] = CoordGen(N,R); PR = THR';
end

phi = log( sum( exp( THR(:,2:end)*theta ) ) );
th = [-phi; theta];
p = exp(THR*th);
%eta = Pr{N}*p;
eta = PR*p;

if nargout > 3
    
rho = zeros(M+1,1); 
for i = 2: 1+M
    %idx = find(THR(i,:));
    idx = logical(THR(i,:));
    r = sum(THR(i,2:1+N));
    
    %rho(i,1) = invTHr{r}(end,:)*log(inv(Pr{r})*eta(idx));
    rho(i,1) = invTHr{r}(end,:)*log(Pr{r}\eta(idx));
    
    %sub_idx{i-1} = idx;
    %sub_eta{i-1} = eta(idx);
    %sub_theta{i-1} = invTHETA{r}*log(inv(P{r})*eta(idx));
    %sub_p{i-1} = exp(THETA{r}*sub_theta{i-1});
end

end
