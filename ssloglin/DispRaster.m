function DispRaster(raw,model)
addpath(genpath('../'));

c = [    0         0    1.0000;
         0    0.5000         0;
    1.0000         0         0;
         0    0.7500    0.7500;
    0.7500         0    0.7500;
    0.7500    0.7500         0;
    0.2500    0.2500    0.2500 ];

N = raw.N;
%T = raw.T;
xs = raw.xs;

if nargin == 2
    n = model.struct.n;
else
    n = length(xs{1});
end

subplot(2,1,1);axis off;

for i = 1:N
    subplot(N,1,i);
    raster(xs{i},c(rem(i-1,7)+1,:));
    
    axis([raw.period(1) raw.period(2) 0 n]);
end

subplot(N,1,1); title(strcat('n=',num2str(n)));
drawnow;
