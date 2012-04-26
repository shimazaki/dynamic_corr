function [BF] = BayesFactor(model,cell_id)

%%%%%% Bayes Factor  
if 0
% Joint Synchrony
clear BF
hoc_idx = find(binary.order == max(binary.order));
delta = zeros(length(hoc_idx),1);
for k = 1: K
	Post0 = mvncdf(delta,fe(hoc_idx,k),fv(hoc_idx,hoc_idx,k));	%Evidence for Positive correlations
	Pred0 = mvncdf(delta,oe(hoc_idx,k),ov(hoc_idx,hoc_idx,k));
    BF(k) = (1 - Post0) / Post0   / ( (1 - Pred0) / Pred0  );
end

stat.BF.filter = BF;
BF(find(isnan(BF))) = 1;
stat.BF.filter_logsum = sum(log2(BF),2);

end

if 0
hoc_idx = find(binary.order == max(binary.order))
for i = 1: length(hoc_idx)
	j = hoc_idx(i);
	for k = 1: K
		switch 1
		case 1 %Positive v.s. Negative 
		delta = 0;
		Post0 = mvncdf(delta,fe(j,k),fv(j,j,k));	%Evidence for Positive correlations
		Pred0 = mvncdf(delta,oe(j,k),ov(j,j,k));
		%BF(i,k) = (1 - Post0) / Post0 / ( (1 - Pred0) / Pred0 );
		
		case 2 %Small v.s. Large
		delta = 2;
		Post0 = 1/2*erf( (delta - fe(j,k))/sqrt(2*fv(j,j,k)) ) ...
				+1/2*erf( (delta + fe(j,k))/sqrt(2*fv(j,j,k)) );
				
		Pred0 = 1/2*erf( (delta - oe(j,k))/sqrt(2*ov(j,j,k)) ) ...
				+1/2*erf( (delta + oe(j,k))/sqrt(2*ov(j,j,k)) );
				
		case 3
		Post0 = mvncdf(delta,fe(j,k),fv(j,j,k)) - mvncdf(-delta,fe(j,k),fv(j,j,k));
		Pred0 = mvncdf(delta,oe(j,k),ov(j,j,k)) - mvncdf(-delta,oe(j,k),ov(j,j,k));
		end
				
		BF(i,k) = (1 - Post0) / Post0   / ( (1 - Pred0) / Pred0  );
        
    end
    
end
stat.BF.filter = BF;

BF(find(isnan(BF))) = 1;
stat.BF.filter_logsum = sum(log2(BF),2);
end
