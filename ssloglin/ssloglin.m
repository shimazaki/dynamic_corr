function param = ssloglin(binary,init)
% This function returns filter/smoother estimate of the spike interactions.
% 
% Shimazaki H., Amari S., Brown E. N., and Gruen S.
% State-space Analysis of Time-varying Higher-order Spike Correlation 
% for Multiple Neural Spike Train Data. 
% PLoS Computational Biology 8(3): e1002385. 
% http://dx.doi.org/10.1371/journal.pcbi.1002385
%
% This program is distributed under the terms of 
% the Creative Commons Attribution License (CC-BY)
% Hideaki Shimazaki, Ph.D.
% http://2000.jukuin.keio.ac.jp/shimazaki

%%%%%%%%%%%%%%% Initialization of Parameters  %%%%%%%%%%%%%%
TRIAL = init.em_max_iteration;

[y,n,TH,miss] = deal(binary.y,binary.n,binary.TH,binary.miss);
[Q,F,G,mu,sig] = deal(init.hyper.Q,init.hyper.F,init.hyper.G,init.hyper.mu,init.hyper.sig);
state = init.state;

[M,K] = size(y);
P = TH';

%u = [zeros(M*2,1), [y(:,1); zeros(M,1)], [y(:,2:K); y(:,1:K-1)]];
%u = [zeros(M*3,1), [y(:,1); zeros(M*2,1)], [y(:,2); y(:,1); zeros(M,1)], ...
%    [y(:,3:K); y(:,2:K-1); y(:,1:K-2)]]; 
u = y(:,1:K);

p = 1; %history
for i = 1: p-1, u = [u; [zeros(M,i), y(:,1:K-i)]]; end
u = [zeros(M*p,1) u];

em_min_logLdiff = init.em_min_logLdiff;
%%%%%%%%%%%%%%%%%%%%%%%  EM-algorithm %%%%%%%%%%%%%%%%%%%%%%%  
trial = 1;
logL_buf = 0; logL_diff = 100;
while (trial < TRIAL+1) && (logL_diff > em_min_logLdiff) 

    %%%% Expectation-step %%%%
    [param.onestep,param.filter,param.smoother] = Estep(F,Q,mu,sig,miss,n,y,P,TH);

    %%%% Maximization-step %%%%
    param.hyper = struct('F',F,'Q',Q,'G',G,'mu',mu,'sig',sig);
    [F,Q,mu,sig] = Mstep(F,Q,mu,sig,param.smoother,state);

    %%%% Marginal Log-likelihood %%%%
    param.stat.logL = MarginalLogLikelihood(param,y,n,TH);
    
    logL = param.stat.logL;
    logL_diff = abs(logL_buf - logL);
    logL_buf = logL;
    
    disp(sprintf('%d logL %g diff %g',[trial logL logL_diff]))
    
    if trial == TRIAL, break; end
    trial = trial + 1;

end % End of Trial

%F, Q, Sig, mu                              
param.stat.dim = M*M*state.F + M*(M+1)/2*state.Q + M*(M+1)/2*0 + M;

param.hyper = struct('Q',Q,'F',F,'G',G,'mu',mu,'sig',sig);     
param.stat = GenStatistics(param,binary);

%%%%%%%%%%%%%%%%%%%  END OF FUNCTION  %%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%  E-step  %%%%%%%%%%%%%%%%%%%%%%%%%%%
function [onestep,filter,smoother] = Estep(F,Q,mu,sig,miss,n,y,P,TH)

[M,K] = size(y);
zi = mu*ones(1,K);

P2 = P(2:end,:);
TH2 = TH(:,2:end);
oneM = ones(1,M);
Ftrans = F';
I = eye(M,M);

%%%%%%%%%%%%%%%%%% Filtering %%%%%%%%%%%%%%%%%%%
oe = zeros(M,K);
ov = zeros(M,M,K);
iov = zeros(M,M,K);

fe = zeros(M,K);
fv = zeros(M,M,K);        

for k = 1: K
    % One-step prediction
    if k == 1
        oe(:,k) = mu;
        ov(:,:,k) = sig;
    else
        oe(:,k) = F*fe(:,k-1);%+G*u(:,k);   %One-step Prediction Expectation
        ov(:,:,k) = F*fv(:,:,k-1)*Ftrans + Q;   %One-step Prediction Variance
    end

    %ov(:,:,k) = (ov(:,:,k)+ov(:,:,k)')/2;
    iov(:,:,k) = inv(ov(:,:,k));
    
    if miss(k) %missing data
    	fe(:,k) = oe(:,k);
    	fv(:,:,k) = ov(:,:,k);
    else
    
    % Posterior mode and variance
    % Newton-Raphson
    z = oe(:,k); %mu;%zi(:,k);
    z_old = z+10;
    count = 1;
    while max( abs(z - z_old) ) > 10^-5 %&& (count < 2) %200
        z_old = z;
        
        % Eta-coordinates
        phi = log( sum( exp( TH2*z ) ) ); 
        %%e = P(2:end,:)*exp(TH*[-phi; z]);
        p = exp(TH*[-phi; z]);
        e = P2*p;
        
        % Fisher Information
        %%FI = P(2:end,:)* (((inv(P(2:end,:))*e)*ones(1,M)) .* P(2:end,:)') - e*e';
        %%FI = P(2:end,:)* ( (exp(TH*[-phi; z])*ones(1,M)) .* (TH*[-e'; eye(M)]) ) ;
        %%FI = P(2:end,:)* ( (EXP_TH*ones(1,M)) .* (TH*[-e'; eye(M)]) ) ;
        FI = P2 * ( (p*oneM) .* TH2 ) - e*e';
        %FI = (FI+FI')/2;

		switch 2
		case 1
        % Newton-Raphson 
        f = z - oe(:,k) - n*ov(:,:,k) * ( y(:,k) - e );
        %f = z - oe(:,k) - G*u(:,k) - n*ov(:,:,k) * ( y(:,k) - e );        
        J = eye(M) + n*ov(:,:,k)*FI;
        z = z - J\f; %inv(J)*f;
        
		case 2
        % Newton-Raphson (faster)
        g = n* ( y(:,k) - e ) - iov(:,:,k) * (z - oe(:,k)) ; 
        H = -n*FI - iov(:,:,k);
        z = z - H\g;

		end
        
        count = count + 1;
        
    end
    if count>199, k, mean( abs(z - z_old) ), disp('ouch'); end
    %count    
    fe(:,k) = z;
    %%fv(:,:,k) = inv( inv( ov(:,:,k) ) + n*FI );
    %%fv(:,:,k) = ( iov(:,:,k) + n*FI ) \ I;  % (faster)
    fv(:,:,k) = inv( iov(:,:,k) + n*FI );   % (fastest)
    
    %fv(:,:,k) = (fv(:,:,k)+fv(:,:,k)')/2;
    
    end
end


%%%%%%%%%%%%%%%%%% Smoothing %%%%%%%%%%%%%%%%%%%
se = zeros(M,K); 
sv = zeros(M,M,K); 

se(:,K) = fe(:,K);      %Initial Value for Smoothing
sv(:,:,K) = fv(:,:,K);     

A = zeros(M,M,K-1);
for h = 1: K-1
    k = K - h;
    %A(:,:,k) = fv(:,:,k)*F'*inv( ov(:,:,k+1) );
    A(:,:,k) = fv(:,:,k)*Ftrans*iov(:,:,k+1);
    
    se(:,k) = fe(:,k) + A(:,:,k)*( se(:,k+1) - oe(:,k+1) );
    sv(:,:,k) = fv(:,:,k) + A(:,:,k)*( sv(:,:,k+1) - ov(:,:,k+1) )*A(:,:,k)'; 
    %sv(:,:,k) = (sv(:,:,k)+sv(:,:,k)')/2;
end

sv2 = zeros(M,M,K-1);
sv2diag = zeros(M,K-1);
for k = 1: K-1
    sv2(:,:,k) = A(:,:,k)*sv(:,:,k+1);
    sv2diag(:,k) = diag(sv2(:,:,k));
end
%sv2 = A(:,:,1:K-1)*sv(:,:,2:K);

fvdiag = zeros(M,K); svdiag = zeros(M,K);
for k = 1: K
    fvdiag(:,k) = diag(fv(:,:,k));
    svdiag(:,k) = diag(sv(:,:,k)); 
    ovdiag(:,k) = diag(ov(:,:,k)); 
end

%sphi = log( sum( exp( TH2*se) ) );
%sp = exp(TH*[-sphi; se]); 
%seta = P2*sp;

onestep = struct('theta',oe,'W',ov,'diagW',ovdiag,'iW',iov);
filter = struct('theta',fe,'W',fv,'diagW',fvdiag);
smoother = struct('theta',se,'W',sv,'diagW',svdiag,'sv2',sv2,'diagsv2',sv2diag,'A',A);
    %'phi',sphi,'eta',seta,'p',sp);

%%%%%%%%%%%%%%%%%%%%%%%  M-step  %%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F,Q,mu,sig] = Mstep(F,Q,mu,sig,smoother,state);

[se,sv,sv2] = deal(smoother.theta,smoother.W,smoother.sv2);
[svdiag,sv2diag] = deal(smoother.diagW,smoother.diagsv2);

[M,K] = size(se);

if (state.Q == 3)       %a full matrix
    v = zeros(M,M);
    for k = 2:K
        v = v + sv(:,:,k)-sv2(:,:,k-1)*F'-F*sv2(:,:,k-1)'+F*sv(:,:,k-1)*F' ...
            + (se(:,k)-F*se(:,k-1))*(se(:,k)-F*se(:,k-1))';
            %+ (se(:,k)-F*se(:,k-1)-G*u(:,k))*(se(:,k)-F*se(:,k-1)-G*u(:,k))';
    end
    Q = v./(K-1);
    Q = (Q+Q')/2;
    
elseif (state.Q == 1)   %a single parameter 120207 works only for F=I.
    v = 0;
    for k = 2:K
        v = v + trace(sv(:,:,k)) + se(:,k)'*se(:,k)...
        -2*trace(sv2(:,:,k-1)) - 2*se(:,k-1)'*se(:,k)...
        +trace(sv(:,:,k-1)) + se(:,k-1)'*se(:,k-1); 
    end
    Q = v/M/(K-1)*eye(M);

elseif (state.Q == 2)   %diagonal
    v = zeros(M,1);
    %v = sum( svdiag(:,2:K)- 2*sv2diag(:,1:K-1)+svdiag(:,1:K-1) ...
    % + se(:,2:K)*se(:,2:K)' -2*se(:,2:K)*se(:,1:K-1)'+se(:,1:K-1)*se(:,1:K-1)', 2 ); 
 
    v = sum( svdiag(:,2:K)- 2*sv2diag(:,1:K-1)+svdiag(:,1:K-1) ...
        + se(:,2:K).^2 -2*se(:,2:K).*se(:,1:K-1)+se(:,1:K-1).^2, 2 ); 
   
    v = v./(K-1);
    Q = diag(v);

elseif (state.Q == 0) %Stationary analysis
	Q = zeros(M,M);
end

if (state.G == 1) && (state.F == 1)
    v1 = 0; v2 = 0;
    w1 = 0; w2 = 0; w3 = 0; w4 = 0;
    for k = 2:K
        v1 = v1 + sv2(:,:,k-1) + se(:,k)*se(:,k-1)';
        v2 = v2 + sv(:,:,k-1) + se(:,k-1)*se(:,k-1)';
    
        w1 = w1 + se(:,k)*u(:,k)';
        w2 = w2 + se(:,k-1)*u(:,k)';
        w3 = w3 + u(:,k)*se(:,k-1)';
        w4 = w4 + u(:,k)*u(:,k)';
    end
    [iv1,iv2,iv3,iv4] = BlockwiseInversion(v2,w2,w3,w4);

    F = v1*iv1 + w1*iv3;
    G = v1*iv2 + w1*iv4;

elseif ((state.F == 1) | (state.F==0.5) )&& (state.G == 0)  
    v1 = 0; v2 = 0;
    for k = 2:K
        v1 = v1 + sv2(:,:,k-1) + se(:,k)*se(:,k-1)';
        v2 = v2 + sv(:,:,k-1) + se(:,k-1)*se(:,k-1)';
    end

    F = v1*inv(v2);
    %G = zeros(M,2*q+1);
    
    if state.F == 0.5
    end

elseif (state.F == 0) && (state.G == 0)
    F = eye(M,M);
    %G = zeros(M,2*q+1);
end

% Initival values update
switch 1
    case 1
        mu = se(:,1);
	case 2
		sig = sv(:,:,1) + (se(:,1) - mu)*(se(:,1) - mu)';
        sig = (sig+sig')/2;

    case 3 %badly behave
        mu = se(:,1);

		sig = sv(:,:,1) + (se(:,1) - mu)*(se(:,1) - mu)';
        sig = (sig+sig')/2;

    case 4 %Equilibrium distribution
        dF = diag(diag(F));
        mu = zeros(M,1);
        vecsig = inv(eye(M*M,M*M)-kron(dF,dF))*reshape(Q,M*M,1);
        sig = reshape(vecsig,M,M)
        
        %sig - (F*sig*F'+Q)
end

%%%%%%%%%%%%%%%%%%%%%%%%%  Q-function  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%function Qfunc = Qfunction


%%%%%%%%%%%%%%%%%%%%%  Marginal Likelihood  %%%%%%%%%%%%%%%%%%%%%%
function logL = MarginalLogLikelihood(param,y,n,TH)

[oe,ov,iov] = deal(param.onestep.theta,param.onestep.W,param.onestep.iW);
[fe,fv] = deal(param.filter.theta,param.filter.W);
[M,K] = size(y);

TH2 = TH(:,2:end);

switch 2
    case 1 %Linear approximation
        logL = 0; 
        for k = 1: K
            phi = log( sum( exp( TH(:,2:M+1)*oe(:,k) ) ) );
            logL = logL + n*( y(:,k)'*oe(:,k) - phi );
        end
    case 2 %Quadratic approximation
        logL = 0;
        for k = 1: K
            phi = log( sum( exp( TH2*fe(:,k) ) ) );
            logL = logL + n*( y(:,k)'*fe(:,k) - phi ) ...
			- 1/2*(fe(:,k)-oe(:,k))'*iov(:,:,k)*(fe(:,k)-oe(:,k)) ...
			 + 1/2*log( det(fv(:,:,k)) ) - 1/2*log( det(ov(:,:,k)) );
        end
    case 3 % Monte Carlo approximation
        logL = 0;
        for k = 1: K
            R = 0; Samples = 500; 
            for i = 1: Samples
                r = mvnrnd(oe(:,k),ov(:,:,k))';
                phi = log( sum( exp( TH(:,2:M+1)*r ) ) );
                R = R + exp( n*(y(:,k)'*r - phi) ) / Samples;
            end
            logL = logL + log(R);
        end
end
      
%%%%%%%%%%%%%%%%%%% Generate Statistics %%%%%%%%%%%%%%%%%%%%%%%%
function stat = GenStatistics(param,binary)
[n K] = deal(binary.n,binary.K);
[dim logL] = deal(param.stat.dim,param.stat.logL);

stat.dim = dim;
stat.logL = logL;

% Information Criteria for Model Selection 
stat.AIC = -2*logL + 2*dim; 
stat.AICc = stat.AIC + 2*dim*(dim+1)/(n-dim-1);
stat.BIC = -2*logL + dim*log(n);

