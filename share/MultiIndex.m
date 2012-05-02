function [subset d] = MultiIndex(N,R,idx)
% This program is distributed under the terms of 
% the Creative Commons Attribution License (CC-BY)
% Hideaki Shimazaki, Ph.D.
% http://2000.jukuin.keio.ac.jp/shimazaki

R = sort(R);

if nargin == 2

	if length(R) == 1
		subset = nchoosek(1:N,R);
		d = length(subset(:,1));

	else
		subset = cell(length(R),1);
		d = 0;
		j = 1;
		for i = R
			subset{j} = nchoosek(1:N,i);
			d = d + length(subset{j}(:,1));
			j = j + 1;
		end
	end


elseif nargin == 3

	ds(1) = 1;
	j = 2;
	for i = 1: length(R) %i = 1: length(R)-1
		ds(j) = nchoosek(N,R(i)); %nchoosek(N,R(i));
		j = j + 1;
	end

	CDS = cumsum(ds);
	rb = max(find(CDS<=idx));		%block that idx belongs to
					
	idxr = 1 + ( idx - CDS(rb) );	%index inside the rth order
									%CDS(r)	%the first index of order r
	ID = nchoosek(1:N,R(rb));
	subset = ID(idxr,:);

end
