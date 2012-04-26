function init = GenInit(init,d)
% This file specifies initial values of the hyper-parameters. 

if ~isfield(init,'em_max_iteration')
	init.em_max_iteration = 200;
end

if ~isfield(init,'state')
    init.state = struct('Q',1,'F',0,'G',0);
end

if ~isfield(init,'em_min_logLdiff')
    init.em_min_logLdiff = 0.1;
end

hyper.mu = zeros(d,1);
hyper.sig = diag( 1*ones(1,d) ); %0.5 for 8 neurons
hyper.F = 1*eye(d);  %0.999
hyper.G = zeros(1,d*10);
hyper.Q = diag(.01*ones(1,d)); % 0.05 in paper; Riehle %0.01^2 for 8 neurons

init.hyper = hyper; 

%Stationary solution
%hyper.mu = zeros(d,1);
%hyper.sig = zeros(d,d);
%vecsig = inv(eye(d*d,d*d)-kron(init.F,init.F))*reshape(init.Q,d*d,1);
%hyper.sig = reshape(vecsig,d,d)

