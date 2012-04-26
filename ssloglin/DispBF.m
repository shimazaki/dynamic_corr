function DispBF(binary,BF,bin)

global ui_slider;
global buf;

buf.binary = binary;
buf.BF = BF;

ui_slider = uicontrol('style','slider','position',[0 0 380,15],...
	'Min',0,'Max',60,'Value',buf.bin,'callback',...
	'buf.bin = get(ui_slider,''value''); clf;buf.bin = round(buf.bin); buf.bin = buf.bin+rem(buf.bin,2); DispBF(buf.binary,buf.BF,buf.bin)' );
	
c = [    0         0    1.0000;
         0    0.5000         0;
    1.0000         0         0;
         0    0.7500    0.7500;
    0.7500         0    0.7500;
    0.7500    0.7500         0;
    0.2500    0.2500    0.2500 ];
    
hoc_idx = find(binary.order == max(binary.order));
K = binary.K;


title(num2str(buf.bin));
subplot(length(hoc_idx)+1,1,1); hold on;
stem(binary.y(hoc_idx,:)','color',[.8 .8 .8],'MarkerSize',0); %axis tight; axis off; 
YLim = get(gca,'YLim'); idx = find(binary.ext); plot(idx,YLim(1)*ones(1,length(idx)),'m^');
axis tight;

for i = 1: length(hoc_idx)
	subplot(length(hoc_idx)+1,1,i+1); hold on; 
	color = c(rem(hoc_idx(i)-1,7)+1,:);
	
	b = [-.5 0 .5 1.0 1.5 2.0];   %Jeffreys
	%b = log10([1 3 20 150]); %Kass&Raftery
	Log10BF = 10.^(b);
	GuideLine = log2(10.^b);

	% Subjective Guidline
	for j = 1: length(GuideLine), 
		line([0, K],[GuideLine(j) GuideLine(j)],'Color',[.8 .8 .8]); 
		%line([0, K],[-GuideLine(j) -GuideLine(j)],'Color',[.8 .8 .8]); 
	end

	%bin = 200; 
 	log2BF_buf = conv(log2(BF(i,:)),ones(1,bin+1));% / (bin+1);
 	log2BF(i,:) = log2BF_buf(bin/2+1:end-bin/2);
 	hold on; bar(log2BF(i,:),'face',[.5,.5,.5],'edgecolor',[.5 .5 .5]); 
	%bar(log2(BF(i,:)),'face',[.5,.5,.5],'edgecolor',[.5 .5 .5]); 
	axis tight;

	
	if i ~= length(hoc_idx), set(gca,'XTickLabel',[]); end
	
	%sigBF = find( log2BF(i,:) >= 1.6610 ); %3.3219
	%subplot(length(hoc_idx)+1,1,1); stem(sigBF,binary.y(hoc_idx(i),sigBF),'Color',color,'MarkerSize',0); 

	subplot(length(hoc_idx)+1,1,1); 
	idx = find(binary.y(hoc_idx(i),:)~=0); 

	for j = 1: length(idx)
		red =  ( min(1,abs(log2(BF(i,idx(j)) ))/GuideLine(3) )) ^2.9 ;
		cr = 1 - [0 red red];
		line([idx(j); idx(j)], length(hoc_idx)-[i i-1],'Color',cr); 
	end
	
	
	if i ~= length(hoc_idx), set(gca,'YTickLabel',[]); end
	YLim = get(gca,'YLim'); idx = find(binary.ext); plot(idx,YLim(1)*ones(1,length(idx)),'m^');
	axis tight;
end



