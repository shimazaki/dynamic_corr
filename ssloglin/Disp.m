function Disp(raw,model)

switch 1
	case 1 % Horizontal
	%w = 380; h = 380; H = 710; o = 20; %for 17inch
    w = 380; h = 380; H = 610; o = 20; %for fisher
	fig(1,1) = figure('Position',[o H w h],'Color','w','Menubar','figure','Toolbar','none','Name','Spike Data');
	fig(1,2) = figure('Position',[1*(w+o)+o H w h],'Color','w','Menubar','figure','Toolbar','none','Name','Joint spike rates');
	fig(1,3) = figure('Position',[2*(w+o)+o H w h],'Color','w','Menubar','figure','Toolbar','none','Name','Interactions');
	fig(1,4) = figure('Position',[3*(w+o)+o H w h],'Color','w','Menubar','figure','Toolbar','none','Name','Hyper-parameter');
	%fig21 = figure('Position',[0 H-h-3*o w h],'Color','w','Menubar','figure','Toolbar','none','Name','Bayes Factor');
    fig(1,5) = figure('Position',[o H-h-4*o 2*w h],'Color','w','Menubar','figure','Toolbar','none','Name','Interactions');
    
	case 2 % Vertical
	w = 400; h = 400; H = 800;
	fig1 = figure('Position',[0 H w h],'Color','w','Menubar','figure','Toolbar','none','Name','Spike Data');
	fig2 = figure('Position',[0 H-h w h],'Color','w','Menubar','figure','Toolbar','none','Name','Posterior');
	fig3 = figure('Position',[0 H-2*h w h],'Color','w','Menubar','figure','Toolbar','none','Name','Hyper-parameter');
	fig4 = figure('Position',[0 H-3*h h 400],'Color','w','Menubar','figure','Toolbar','none','Name','Bayes Factor');

end

figure(fig(1,1)); DispRaster(raw,model); 
figure(fig(1,2)); DispBinary(raw,model); 
figure(fig(1,3)); %DispEstimates(model,'filter'); 
			  	  DispEstimates(model,'smoother');
figure(fig(1,4)); DispHyperparam(model.param.hyper);
figure(fig(1,5)); DispConnections(model,2);

global buf;
%figure(fig4); DispBF(model.binary,model.param.stat.BF.filter,buf.bin);

drawnow;
