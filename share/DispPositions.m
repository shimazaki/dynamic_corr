function DispPositions(position,cell_id,c)
% Display Cell Positions

[N,a] = size(position);

if nargin <= 1, cell_id = 1:N; end

if nargin <= 2, c = [0 0 0]; end



hold on;
plot(position(1:N,1),position(1:N,2),'.','Color',c,'MarkerSize',10); 
%plot(xr,yr,'ko','MarkerSize',5); 

for i = 1: N
    text( position(i,1)+4,position(i,2),num2str( cell_id(i) ) );
end

set(gca,'TickDir','out');
