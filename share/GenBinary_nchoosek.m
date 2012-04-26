function binary = GenBinary(raw,struct)
%

[N,R,D,n] = deal(struct.N,struct.R,struct.D,struct.n);
[xs] = deal(raw.xs);
K = ceil(raw.T / D);

ts = raw.period(1); tf = raw.period(2);
dt = raw.D; 


pattern = zeros(K,N,n);
for trial = 1: n
	for i = 1: N
		xsb = xs{i}{trial};
		xsb = xsb(find( (ts < xsb ) .* (xsb < tf) )) - ts;

		pattern(ceil(( xsb )/D ),N-i+1,trial) = 1; %original
		%pattern(ceil(( xs{i}{trial} )/D),N-i+1,trial) = 1;
	end
end

Rp = R(R > 0);
Rn = abs(R(R < 0));
[idx_buf dp] = MultiIndex(N,Rp);
[idx_buf dn] = MultiIndex(N,Rn);

d = dp + dn;
y = zeros(d,K);

for trial = 1: n
	for k = 1: K
		subset = find(pattern(k,:,trial));
		if isempty(subset) %sum(pattern(k,:,trial)*pattern(k,:,trial)') == 0
		else			
			[idx idxs] = invMultiIndex(N,Rp,subset);
			y(idxs,k) = y(idxs,k) + 1/n;
		end

		subset = find(~pattern(k,:,trial));
		if isempty(subset)
		else
        	[idx idxs] = invMultiIndex(N,Rn,subset); %need to be double checked!!
        	y(dp+idxs,k) = y(idxs,k) + 1/n;
		end
	end	
end
dp,dn

%%%%%%%%%%%%%%%%%%%%% binary %%%%%%%%%%%%%%%%%%%%%%%
binary.n = n;
binary.d = length(y(:,1));
binary.K = K;

binary.y = y;

binary.y_mean = mean(y,2);

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


%
N,R
TH = CoordGen_nchoosek(N,R); P = TH';
binary.TH = TH;

