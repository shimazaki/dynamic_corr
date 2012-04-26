function DispConnections(model,r)
% DispConnections displays snapshots of spike correlations
%

theta = model.param.smoother.theta;
N = model.struct.N;
K = model.binary.K;

%%%%%%%%%%%%%%%%%%%% Initial Settings %%%%%%%%%%%%%%%%%%%%
R = 2;
SHOTS = 5;

for r = 1 : N
    [THr{r}]= CoordGen(r);
	Pr{r} = THr{r}';
	invTHr{r}=inv(THr{r});
end

for k = 1: K
    [thetah(:,k),etah(:,k),ph(:,k)] ...
        = CoordTrans(theta(:,k),N,R,Pr,THr,invTHr);
end

hoc_idx = find(model.binary.order == max(model.binary.order))';
IDX = model.binary.id;

rate = etah(2:N+1,:);
rho = theta;
rho(isinf(rho)) = 100; 

%%
[M,K] = size(rho);
rad = pi/2 - ([0:N-1]+.0)* 2*pi / N;
xr = cos(rad);
yr = sin(rad);


ts = round(linspace(100,K,SHOTS)); 

maxr = max(max(rate(1:N,:)))
minr = min(min(rate(1:N,:)))
            
maxb = max(max(rho(hoc_idx,:)));
minb = -min(min(rho(hoc_idx,:)));             

for shot = 1: length(ts)
time = ts(shot); 
	lw = 0.8;
	if 1
		set(gca,'drawmode','fast');
		subplot(2,ceil(length(ts)/2),shot);
	else
		FRAME(shot) = getframe;
    	end

	for i = N+1: M
		idxes = find(IDX(i,:));
		idx1 = idxes(1);
		idx2 = idxes(2);
		rho_posi = rho(i,time)*(sign(rho(i,time))+1)/2;
		rho_nega = rho(i,time)*(sign(rho(i,time))-1)/2;
                
        if max(abs(rho(i,:))) ~= 0
                if rho_nega ~= 0
                    c = 1-1*rho_nega/minb;
                    patch([xr(idx1),xr(idx2)],[yr(idx1),yr(idx2)],[0 1 0],...
                        'EdgeColor',[c,c,1],'Linewidth',lw*rho_nega);
                end
                if rho_posi ~= 0
                    c = 1-1*rho_posi/maxb;
                    patch([xr(idx1),xr(idx2)],[yr(idx1),yr(idx2)],[0 1 0],...
                        'EdgeColor',[1,c,c],'Linewidth',lw*rho_posi);
                end

        end
    	end
          
	% Firing Rate
 	for i = 1: N
        c = (rate(i,time) - minr) / (maxr - minr)+0.001;
        c_idx = ceil(c*63);
        cmap = flipud(hot);
        hold on; plot(xr(i),yr(i),'yo','MarkerSize',5,...
            'MarkerEdgeColor',cmap(c_idx,:),'MarkerFaceColor',cmap(c_idx,:))
        
        if shot == 1
            text(1.1*xr(i),1.1*yr(i),num2str(i)); 
        end
    end
    cmap
    %colorbar(cmap)
	axis square; axis off; 

    title(num2str(time));
    drawnow; 
    

end
