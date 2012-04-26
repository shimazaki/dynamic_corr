function [cmap] = DispInteractions(N,R,theta,r,position,th_max)

[TH IDX order IA] = CoordGen(N,R);
%[TH IDX order] = CoordGen(N,R);

d = length(theta);

x = position;

idx_rs = find(sum(IA,2) == r);
%idx_rs = find(sum(IDX(2:end,:),2) == r);

%th_max = 0;
%for i = 1: length(idx_rs)
%	idx_r = idx_rs(i);
%	th = theta(idx_r);
%	th_max = max(th_max,abs(th));
%end

for i = 1: N
    %hold on; plot(x(i,1),x(i,2),'yo','MarkerSize',3,...
    %       'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1])
end

for i = 1: length(idx_rs)
	idx_r = idx_rs(i);
	cell_id = find(IA(idx_r,:));
    %cell_id = find(IDX(idx_r,:));

	x_mean = mean(x(cell_id,:),1);

	th = theta(idx_r);
	%c = 0.9;
	%c = min(1,abs(th)/1);
	c = min(1, abs(th) / th_max);

	lw = 2;0.5;1.2;
	if th >= 0
		linecolor = [1,1-c,1-c];
		linewidth = lw;%lw*abs(th);
	else
		linecolor = [1-c,1-c,1];
		linewidth = lw;%lw*abs(th);
	end
	
	for j = 1: length(cell_id)
		cell = cell_id(j);
		x_j = x(cell,:);

		if abs(th) < 0.2
		else
		hold on; patch([x_mean(1),x_j(1)],[x_mean(2),x_j(2)],[0 1 0],...
 			'EdgeColor',linecolor,'Linewidth',linewidth);
		end
	end 
end


cs = linspace(0,1,32);
for i = 1: 32
	cmap(33-i,:) = [1-cs(i),1-cs(i),1];
end
for i = 33: 64
	cmap(i,:) = [1,1-cs(i-32),1-cs(i-32)];
end

colormap(cmap); %colorbar; set(gca,'CLim',[-th_max th_max]);
