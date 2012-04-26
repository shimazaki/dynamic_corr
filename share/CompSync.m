function sync = CompSync(raw,D,r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Synchronous spikes

N = raw.N;
xs = raw.xs;
n = length(xs{1});

%D = 0.05;
%r = 3;

for i = 1: n
    for j = 1: N
        xsr{j}{i} = ceil( xs{j}{i} /D ) * D;
    end
end

[subset] = MultiIndex(N,r);
[ROW,COL] = size(subset);

for i = 1: n  
    for j = 1: ROW
        xsb{j}{i} = xsr{subset(j,1)}{i};
        for k = 1: COL
            xsb{j}{i} = intersect(xsb{j}{i}, xsr{subset(j,k)}{i});
        end
    end
end


for i = 1: n
	sync{i} = [];
	for j = 1: ROW
		sync{i} = union(sync{i}, xsb{j}{i});
	end
end
