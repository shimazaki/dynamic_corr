function xs = GenSpike(p,n,dt)
%clear all;
%function theta = SignalGen(T,dt,Rgen) 
% xs{i} 	cell array, superimposed spikes

RandStream.setDefaultStream ...
     (RandStream('mt19937ar','seed',sum(1000*clock)));
%rand('state',sum(100*clock))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter Settings
[Md,K] = size(p);
N = log2(Md);
M = Md-1; 

%[THgen IDX] = CoordGen(N,1);    %Coordtrans Matrix original

%%%%%%%%%%%%%%%  Generation of Spikes  %%%%%%%%%%%%%%%%%%
X = sparse(M+1,K);
Xs = cell(1,n);
xss = cell(N,n); 

for i = 1: n
    Xs{i} = sparse(M+1,K); 
    for k = 1: K
        c = cumsum(p(:,k));

        idx = 1 + sum(c<rand);
        X(idx,k) = X(idx,k) + 1;
        Xs{i}(idx,k) = Xs{i}(idx,k) + 1;
        
        %idx2 = find(IDX(idx,:)); % original

		%added 110927
		if idx == 1
			idx2 = [];
		else
			idx2 = MultiIndex(N,1:N,idx-1);
		end
		%addition ends here

        for j = 1: length(idx2)
            xss{idx2(j),i}= [[xss{idx2(j),i}], k*dt];
        end
    end
    
end

xs = cell(1,N);
for i = 1: N, xs{i} = xss(i,1:n); end


%clear xss

%X = zeros(2^N,K);
%for i = 1: n
%    X = X + Xs{i}(:,:);
%end


