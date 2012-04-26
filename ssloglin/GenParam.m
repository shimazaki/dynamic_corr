function [gen,rho] = GenParam(K,N,Rgen)
%gen.
% theta, eta, p
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter Settings
M = 2^N-1; 
dt = 0.001;
T = K*dt;
    
Mgen = 0; 
for r = 1: Rgen
	Mgen = Mgen + nchoosek(N,r);
end
        
theta = zeros(Mgen+1,K);
for i = 2:N+1
    theta(i,:) = StochProcess_Gauss(T,dt,-15/N,3,10+5*(rand-.5));
end
for i= N+2:Mgen+1
    %theta(i,:) = StochProcess_Gauss(T,dt,5.0/N*(0.5*rand+1),1.3+0.2*(rand-.5),8);
    theta(i,:) = StochProcess_Gauss(T,dt,0.5/N*(0.5*rand+1),2,10);
end

j = 2;
for i = 1: N
    theta(i,:) = StochProcess_Gauss(T,dt,-30/N,3,60);
    j = j + 1;
end

for r = 2: N
	for i = 1: nchoosek(N,r)
		theta(j,:) = StochProcess_Gauss(T,dt,6/nchoosek(N,r),3,50);
		j = j + 1;
	end
end
                
        
[THgen IDX] = CoordGen(N,Rgen); Pgen = THgen';    %Coordtrans Matrix 
%load CoordTransMatrix10 Pr THr invTHr   %for up to 10th-order correlations

gen = struct;

eta = zeros(Mgen+1,K); p = zeros(Mgen+1,K); rho = zeros(Mgen+1,K);
if nargout < 2
    for k = 1: K
        %[theta(:,k),eta(:,k),p(:,k),rho(:,k)] = CoordTrans(theta(2:end,k),N,Rgen,Pr,THr,invTHr);
		[theta(:,k),eta(:,k),p(:,k)] = CoordTrans(theta(2:end,k),N,Rgen);
    end
    
    gen.theta = theta;
    gen.eta = eta;
    gen.p = p;
else    
    for k = 1: K
        [theta(:,k),eta(:,k),p(:,k)] = CoordTrans(theta(2:end,k),N,Rgen,Pr,THr,invTHr);
    end
end


