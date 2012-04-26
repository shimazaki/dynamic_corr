function [TH, IDX, order TABLE] = CoordGen(N,RS)
% TH is the matrix that converts the theta-coordinate into the p-coordinate.
% IDX contains all configurations of the binary state.
% p-coordinate: joint probability distribution. 
% thteta-coordinate: log-linear parameters.
% eta-coordinate: expectation of the joint spikes.
% P is the matrix that converts the p-coordinate into the eta-coordinate.
% P is obtained as TH'.

% 2011/6/6 R -> RS, now the function accepts column vetor of orders.


RSp = RS(RS > 0);
RSn = RS(RS < 0);

IDX = [];

if isnan(RSp) == 0
    [THp IDXp orderp] = CoordGen_Subfunction(N,RSp);
	IDX = [IDXp];
else
    THp = []; IDXp = []; orderp = [];
end

if isnan(RSn) == 0
    [THn IDXn ordern] = CoordGen_Subfunction(N,abs(RSn));
    THn = flipud(THn);
    IDXn = 1-flipud(IDXn); 
	
	IDX = [IDX; IDXn];
else
    THn = []; IDXn = []; ordern = [];
end

TH = [ones(2^N,1) THp(:,2:end) THn(:,2:end)];
order = [orderp ordern];% ordern


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [TH IDX order] = CoordGen_Subfunction(N,RS)
% fast method to create indexes.
str =  dec2bin(0:2^N-1,N)';
IDX = sscanf(str,'%1d',[N,2^N])';

% Sorting in the order of the number of joint spikes
c = sum(IDX,2);
[buf,order1] = sort(c);
IDX = IDX(order1,:);
IDX = fliplr(IDX);

[buf order] = sort(order1);

%%% Selction up to Rth order
% Number of elements

RS = sort(RS);

% for using only selected orders
SelectedRowIDXs = [];
M = 0;
MSelected = 0;
for r = 1: N
    for i = 1:length(RS)
        if RS(i) == r
            SelectedRowIDXs =[SelectedRowIDXs M+1:M+nchoosek(N,RS(i))];
            MSelected = MSelected + nchoosek(N,RS(i));
        end
    end
    M = M + nchoosek(N,r);
end
SelectedRowIDXs = SelectedRowIDXs + 1;

% Obtain the matrix TH by using the sate IDX. 
TH = ones(2^N,MSelected+1);
for i = 1: length(SelectedRowIDXs)               %Compute each column of TH
    SI = SelectedRowIDXs(i);                    %selected index
	idx = find(IDX(SI,:));                  %active binary in ith row
	for j = 1: length(idx)
		%column vector that active neurons in idx(j) join
        TH(:,i+1) = and( TH(:,i+1), IDX(:,idx(j)) ) ;
    end
end
