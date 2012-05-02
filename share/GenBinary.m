function binary = GenBinary(raw,struct)
% This program is distributed under the terms of 
% the Creative Commons Attribution License (CC-BY)
% Hideaki Shimazaki, Ph.D.
% http://2000.jukuin.keio.ac.jp/shimazaki

[N,R,D,n] = deal(struct.N,struct.R,struct.D,struct.n);
[xs] = deal(raw.xs);

[TH IDX order] = CoordGen(N,R); P = TH';

ts = raw.period(1); tf = raw.period(2);
K = ceil( (tf-ts) / D);
dt = raw.D; %SHIFT = D/dt;

for shift = 1: 1 %D/dt %shift-average

pattern = zeros(K,N,n);
for trial = 1: n
	for i = 1: N
		xsb = xs{i}{trial}+(shift-1)*dt;
		xsb = xsb(find( (ts < xsb ) .* (xsb < tf) )) - ts;

		pattern(ceil(( xsb )/D ),N-i+1,trial) = 1; %original
		%pattern(ceil(( xs{i}{trial} )/D),N-i+1,trial) = 1;
	end
end

X = sparse(2^N,K); 
base = 2.^fliplr(0:N-1)';
for trial = 1: n
	%idx = bin2dec( num2str(pattern(:,:,trial)) ) + 1;
	idx = pattern(:,:,trial) * base + 1; 	%bin to dec
	idx = order(idx);
	
	%for k = 1: K, X_full(idx(k),k) = X_full(idx(k),k) + 1; end
	X = X + sparse(idx',[1:K],1,2^N,K); %faster way to obtain X
end

y_buf(:,:,shift) = full((P*X)/n);

end

y = mean(y_buf,3); %y = full((P*X)/n);
y = y(2:end,:);

%%%%%%%%%%%%%%%%%%%%% binary %%%%%%%%%%%%%%%%%%%%%%%
%%binary = struct;
binary.n = n;
binary.d = length(y(:,1));
binary.K = K;

binary.y = y;

binary.y_mean = mean(y,2);
binary.p_mean = full(sum(X,2)/K/n);

% the next two variables has 2^N elements.
binary.id = full(IDX(2:end,:));
binary.order = full( sum(IDX(2:end,:),2) );

%binary.y_order = full( sum(IDX(2:end,:),2) );

ext = zeros(1,K); 
if isfield(raw,'ext')	
	idx = ceil( (raw.ext - ts + dt/2) / D);
	idx = idx(idx>0);
	ext( idx ) = 1;	
end 
binary.ext = sparse(ext);

miss = zeros(1,K); 
if isfield(raw,'miss')
	miss(ceil( (raw.miss - ts + dt/2) / D)) = 1;
end
binary.miss = sparse(miss);

binary.X = X;
binary.TH = TH;

