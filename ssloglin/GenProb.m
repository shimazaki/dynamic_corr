function [p eta theta] = GenProb
path(path,'~/lab/tool/');
T = 0.5; dt = 0.001; mu = -3; s = 0.5; g = 50;
L = T/dt; tt = [dt:dt:T];

switch 3
	case 1
		y(1,:) = ones(1,L);
		y(2,:) = StochProcess_Gauss(T,dt,-3,s,g);
		y(3,:) = StochProcess_Gauss(T,dt,-3,s,g);
		y(4,:) = StochProcess_Gauss(T,dt,1,1,100);
	case 2 %Rate-covariate No-sync
        w= 0.01;
		y(1,:) = ones(1,L);
		y(2,:) = 2*w^2./(w^2+(tt-T/2).^2)-3;
		y(3,:) = 2*w^2./(w^2+(tt-T/2).^2)-3.5;
		y(4,:) = 0*ones(1,L);
	case 3 %Rate-cosntant Sync-nonstationary
        w = 0.01;
        y(1,:) = ones(1,L);
		y(2,:) = -3*ones(1,L); StochProcess_Gauss(T,dt,mu,0,g);
		y(3,:) = -3.5*ones(1,L);StochProcess_Gauss(T,dt,mu,0,g);
		y(4,:) = -5*w^2./(w^2+(tt-T/2).^2);
end

% Rectification
y(1,:) = y(1,:).*(y(1,:)<0);
y(2,:) = y(2,:).*(y(2,:)<0);
rho = y;

% Backward calculation
p = zeros(4,L); eta = p; theta = p;
for t = 1: L
	disp(t);
	[p(:,t), eta(:,t) ,theta(:,t)] = CoordTrans2(y(:,t));
end



function [p,eta,theta] = CoordTrans2(rho)

N = round(log2(length(rho)));

eta = zeros(2^N,1);
p = zeros(2^N,1);
theta = zeros(2^N,1);

eta(1) = exp(rho(1));

P{1} = [1 1;
        0 1];
P{2} = [1 1 1 1;
        0 1 0 1;
        0 0 1 1;
        0 0 0 1];
r = 1;
[eta(2) RHO(2,:)]= NewtonRaphson(P{r},eta(1),rho(2));
[eta(3) RHO(3,:)]= NewtonRaphson(P{r},eta(1),rho(3));

r = 2;
[eta(4) RHO(4,:)]= NewtonRaphson(P{r},eta(1:3),rho(4));

p = inv(P{N})*eta;
theta = inv(P{N}')*log(inv(P{N})*eta);


function [eta_l RHO]= NewtonRaphson(P, eta_in, rho_l)
invP = inv(P);
invPp = inv(P');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial value for a Newton-Raphson Method
ND = 10^3;
xs = linspace(10^(-10),1-10^(-10),ND);  %Candidate eta_l
ys = zeros(1,ND);
for i = 1: ND
    ys(i) = theta_l(P,invP,invPp,[eta_in; xs(i)]);
end
%subplot(3,6,7); plot(xs,ys); hold on; plot(xs,0,'k-'); axis fill; hold off; drawnow;
xs = xs(~isnan(ys));
ys = ys(~isnan(ys));
idx = find(abs(diff( ys -rho_l < 0 )));
if isempty(idx) == 1
    x = NaN;
else
    x = xs(idx);
end
RHO = [min(ys), max(ys)];
eta_l = x;

if isnan(x) ~= 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newton-Raphson 
eta = [eta_in; x];
dx = 10^-5;
for i = 1: 500
    if isnan( theta_l(P,invP,invPp,[eta_in; x+dx] ) ) ~= 1
        eta(end) = x;
        df = ( theta_l(P,invP,invPp,[eta_in; x+dx]) - theta_l(P,invP,invPp,eta) )/dx;
        x = x - (theta_l(P,invP,invPp,eta)-rho_l)/df;
    end
end
eta_l = eta(end);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Return the highest theta
function y = theta_l(P,invP,invPp,eta)
p = invP*eta;
%if sum(p <= 0) + sum(p > 1) > 0
%    y = NaN;
%else
%    y = invPp(end,:)*log(p);
%end

flag = 0;
for i = 1: length(p)
    if p(i) <= 0 
        flag = 1;
    elseif p(i) > 1
        flag = 1;
    end
end

if flag == 1
    y = NaN;
else
    y = invPp(end,:)*log(p);
end


