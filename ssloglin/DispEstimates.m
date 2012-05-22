function DispEstimates(model,type,gen)

c = [    0         0    1.0000;
         0    0.5000         0;
    1.0000         0         0;
         0    0.7500    0.7500;
    0.7500         0    0.7500;
    0.7500    0.7500         0;
    0.2500    0.2500    0.2500 ];
    
cf = 1.95996; %95
%cf = 2.57583; %99 %  3; 2.32635; 2.80703; 

mode = model.param.(type).theta; 
interval = model.param.(type).diagW;

N = model.struct.N;
R = model.struct.R;

K = model.binary.K;
d = model.binary.d;

%order = model.binary.order; 
ext = model.binary.ext; 

%%%%%%%%%%%%%%%% Credible Interval %%%%%%%%%%%%%%
if 1
    
j = 1;
for i = 1: 0;%d
	order = length(MultiIndex(N,R,i));
	subplot(length(R), 1, find(R==order)); 
    
	hold on; line([1:K; 1: K],...
	        [mode(j,:)+cf*sqrt(interval(j,:)); mode(j,:)-cf*sqrt(interval(j,:))]...
            	,'Color',[7 7 7]/8,'LineWidth',2 );       
    j = j + 1;
end

j = 1;
for i = 1: d
	order = length(MultiIndex(N,R,i));
	subplot(length(R), 1, find(R==order)); 

	hold on; line([1 K],[0 0],'Color',[.6 .6 .6]);

    hold on; plot(1:K,mode(j,:)+cf*sqrt(interval(j,:)),'Color',c(rem(i-1,7)+1,:));
    hold on; plot(1:K,mode(j,:)-cf*sqrt(interval(j,:)),'Color',c(rem(i-1,7)+1,:));
    axis tight;     
    
    j = j + 1;
end

end

%%%%%%%%%%%%%%%%%%%% Mode %%%%%%%%%%%%%%%%%%%%%%%
j = 1;
for i = 1: d
	order = length(MultiIndex(N,R,i));
	subplot(length(R), 1, find(R==order)); 

    hold on;
    if nargin > 2
        plot(gen.param.theta(j+1,:)','--','Color',c(rem(i-1,7)+1,:)); % Underlying 
    end
    
	plot(mode(j,:),'Color',c(rem(i-1,7)+1,:));
    
    title( strcat('r=',num2str(order)) );
    
	j = j + 1; 
	
	YLim = get(gca,'YLim'); 
	idx = find(ext); plot(idx,YLim(1)*ones(1,length(idx)),'m^');
	axis tight;
	if i ~= d
		set(gca,'XTickLabel',[]);
	end
end



