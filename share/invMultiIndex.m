function [idx idxs] = invMultiIndex_nchoosek(N,R,subset)
% This program is distributed under the terms of 
% the Creative Commons Attribution License (CC-BY)
% Hideaki Shimazaki, Ph.D.
% http://2000.jukuin.keio.ac.jp/shimazaki

%nCk =  @(n,k) prod(n-(k-1):n)/prod(1:k);
R = sort(R);

r = length(subset);

idxR = 0;
for j = R(R < r)
	idxR = idxR + nCk(N,j);%nchoosek(N,i);
end

idx = idxR + invMultiIndex_Single(N,R,subset);

if nargout == 2
idxs = [];

k = length(subset);							%the k-subset
% recursively create i-subsets of k elements (i <= k). 
% i-subsets are specified by interaction order R (not larger than k).		
R_not_larger_than_k = R( R <= k );

m = 1;
for i = R_not_larger_than_k
	isub = nchoosek(subset,i);				% i-subset of id (an element of k-subset). 

	r = length(isub(1,:));
	
	idxR = 0;
	for j = R(R < r)
		idxR = idxR + nCk(N,j); %nchoosek(N,i);
	end

	for l = 1: length(isub(:,1))									 
		idxs(m) = idxR+invMultiIndex_Single(N,R,isub(l,:)); % index of the subset `isub' in (N,R) model.		
		m = m + 1;
	end

end

end


function idxr = invMultiIndex_Single(N,R,subset)
%nchoosek(n,k) is very slow.
%nCk =  @(n,k) prod(1:n)/(prod(1:n-k)*prod(1:k));
%nCk =  @(n,k) prod(n-(k-1):n)/prod(1:k);

r = length(subset);

subset = [0  subset];
idxr = 0;

for i = 1 : r
	ii = i + 1;
	for j = subset(ii-1)+1 : subset(ii)-1
		idxr = idxr + nCk(N-j,r-i); %nchoosek(N-j,r-i);
		%[i j N-j r-i nchoosek(N-j,r-i)] %sanity check
	end
end

idxr = idxr + 1;

%idx = idx + idxr;

function y = nCk(n,k)
y = prod(n-(k-1):n)/prod(1:k);
