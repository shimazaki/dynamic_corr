function DispRaster(raw,model)

c = [    0         0    1.0000;
         0    0.5000         0;
    1.0000         0         0;
         0    0.7500    0.7500;
    0.7500         0    0.7500;
    0.7500    0.7500         0;
    0.2500    0.2500    0.2500 ];

N = model.struct.N;
R = model.struct.R;

K = model.binary.K;
%order = model.binary.order;
d = model.binary.d;

j = 1; 
for i = 1: d
	order = length(MultiIndex(N,R,i));
	subplot(length(R), 1, find(R==order)); 
    
    title( strcat('r=',num2str(order)) );

	hold on; stairs(1:K,model.binary.y(j,:),'Color',c(rem(j-1,7)+1,:));
	axis tight;
  	j = j + 1;
end

drawnow;
