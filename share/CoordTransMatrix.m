function TH = CoordTransMatrix(N,R)

IDs = MultiIndex(N,R);

m = 1;

for k = 1: N

	ids = MultiIndex(N,k);							% k-subset of N elements

	for j = 1: length(ids(:,1))			

		ksub = ids(j,:);							% one sample of the k-subset

		% recursively create i-subsets of k elements (i <= k). 
		% i-subsets are specified by interaction order R (not larger than k).		
		R_not_larger_than_k = R( R <= k );	
	
		for i = R_not_larger_than_k

			isub = nchoosek(ksub,i);				% i-subset of id (an element of k-subset). 

			for l = 1: length(isub(:,1))
										 
				idx = invMultiIndex(N,R,isub(l,:)); % index of the subset `isub' in (N,R) model.

				TH(m, idx) = 1;

			end			

		end
		
		m = m + 1;
	end

end
