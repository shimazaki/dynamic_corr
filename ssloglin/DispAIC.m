function DispABIC(model)


R = size(model(:));

for i = 1: R
	ABIC(i) = model(i).param.stat.ABIC;
end

plot(1:R,ABIC,'o-','MarkerSize',8);
