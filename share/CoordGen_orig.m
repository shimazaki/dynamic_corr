function [TH, IDX, order] = CoordGen(N,R)
% TH is the matrix that converts the theta-coordinate into the p-coordinate.
% IDX contains all configurations of the binary state.
% p-coordinate: joint probability distribution. 
% thteta-coordinate: log-linear parameters.
% eta-coordinate: expectation of the joint spikes.
% P is the matrix that converts the p-coordinate into the eta-coordinate.
% P is obtained as TH'.

if (nargin == 1)
    R = N;
end


IDX = sparse(2^N,N);

% slow
%clear i j
%for i = 1 : 2^N
%	str = dec2bin(i-1,N);               %starting from 00..0
%	for j = 1: N, IDX(i,j) = str2double(str(j)); end
%    IDX(i,:) = str2num(sprintf('%c,',str));
%    IDX(i,:) = sscanf(str,'%1d',[1,N]);
%end

% faster method to creat indexes.
%for i = 1 : 2^N
%    str(:,i) = dec2bin(i-1,N);
%end

str =  dec2bin(0:2^N-1,N)';
IDX = sscanf(str,'%1d',[N,2^N])';

% Sorting in the order of the number of joint spikes
c = sum(IDX,2);
[buf,order1] = sort(c); 
IDX = IDX(order1,:);
IDX = fliplr(IDX);

[buf order] = sort(order1);

M = 0; 
for r = 1: R
    M = M + nchoosek(N,r);
end
 
% Obtain the matrix P by using the sate IDX. 
TH = sparse(2^N,M+1); 
TH = ones(2^N,M+1);
for i = 2: M+1                              %Compute each row of P
	idx = find(IDX(i,:));                   %active binary in ith row
    %%TH(:,i) = ones(2^N,1);
	for j = 1: length(idx)
		%vector that neuron idx(j) joins
        TH(:,i) = and( TH(:,i), IDX(:,idx(j)) ) ;
    end
end

%P = TH';

% Obtain the matrix P by using the sate IDX. 
%P = sparse(M+1,2^N); 
%P(1,:) = ones(1,2^N);
%for i = 2: M+1                              %Compute each row of P
%	idx = find(IDX(i,:));                   %active binary in ith row
%    P(i,:) = ones(1,2^N);
%	for j = 1: length(idx)
%		%P(i,:) = P(i,:) .* IDX(:,idx(j))' ;   %vector that neuron idx(j) joins
%        P(i,:) = and(logical(P(i,:)), IDX(:,idx(j))') ;
%    end
%end


if nargout > 3
    [buf,b2p]=sort(order1);     
    
    %Table for Fisher Information Matrix
    TABLE = zeros(M,M);
    for i = 1: M %2^N
        for j = 1: M
            idx = or(IDX(i+1,:),IDX(j+1,:));
            idx2 = bin2dec(num2str(idx));

            TABLE(i,j) = b2p(idx2+1) - 1;
            %TABLE(i,j) = idx2 - 1;
        end
    end
    %TABLE = TABLE(2:end,2:end);
end
