function [TH] = CoordGen(N,R)
% TH is the matrix that converts the theta-coordinate into the p-coordinate.
% IDX contains all configurations of the binary state.
% p-coordinate: joint probability distribution. 
% thteta-coordinate: log-linear parameters.
% eta-coordinate: expectation of the joint spikes.
% P is the matrix that converts the p-coordinate into the eta-coordinate.
% P is obtained as TH'.

% 2011/6/6 R -> RS, now the function accepts column vetor of orders.


if min(R) > 0
    
    TH = CoordTransMatrix(N,R);
    [m,n] = size(TH);
    TH = [ones(2^N,1) [zeros(1,n); TH] ];

else
    Rp = R(R > 0);
    Rn = R(R < 0);
	TH = [];
    
	if isnan(Rp) == 0
		[THp] = CoordTransMatrix(N,Rp);
		[mp,np] = size(THp);
		TH = [TH [zeros(1,np); THp]];
    end

	if isnan(Rn) == 0
		[THn] = CoordTransMatrix(N,abs(Rn));
		[mn,nn] = size(THn);
		TH = [TH flipud([zeros(1,nn); THn])];
	end

    TH = [ones(2^N,1)  TH];
end
