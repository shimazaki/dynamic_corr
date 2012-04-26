function DispConnections(model,r)
% DispConnections displays snapshots of spike correlations
%

theta = model.param.smoother.theta;
N = model.struct.N;
R = model.struct.R;
K = model.binary.K;

%%%%%%%%%%%%%%%%%%%% Initial Settings %%%%%%%%%%%%%%%%%%%%
SHOTS = 5;

rad = pi/2 - ([0:N-1]'+.0)* 2*pi / N;
position = [cos(rad) sin(rad)];

ks = floor(linspace(1,K,SHOTS));

rs = setdiff(R,1);
th_max = max(max(abs( theta(:,ks) ) ));

for i = 1:length(ks)
    k = ks(i);
    %[thetah,etah,ph] ...
    %    = CoordTrans(theta(:,k),N,R,Pr,THr,invTHr);
    
    for j = 1: length(rs)
    r = rs(j);
    subplot(length(rs),SHOTS,(j-1)*SHOTS+i)
    DispInteractions(N,R,theta(:,k),r,position,th_max);
    title(strcat('r=',num2str( r ),', ',num2str(k),'[\Delta]'));
    axis square; axis off; box off;    
    end
end
