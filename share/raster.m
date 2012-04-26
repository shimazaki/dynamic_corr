function raster(x,c,Properties)
% 
% raster(xs)
% xs cell array
% c color
% for population raster([xs{:}])

if nargin == 1
    c = 'k';
end

if iscell(x) == 1
    xs = x;
    [m,n] = size(x);
else
    n = 1;
    xs{1} = x;
end

hold on;
for i = 1: n
    len = length(xs{i}); 
    p = line([xs{i}; xs{i}], [(n-i)*ones(1,len); (n-i+1)*ones(1,len)]); %n=1 top
    set(p,'Color',c);
    %eval( strcat('set(p',Properties,')') )
end

set(gca,'TickDir','out');
%set(gca,'TickLength',[0 0]);
%axis([Ts Tf 0 n]);
