function GenSyncSpikes(xs,R)



for i = 1: n
	sync12{i} = intersect(ceil((xsb{1}{i}-ts)/D), ceil((xsb{2}{i}-ts)/D) ) *D;
	sync123{i} =intersect(ceil((xsb{1}{i}-ts)/D), ...
			intersect(ceil((xsb{2}{i}-ts)/D), ceil((xsb{3}{i}-ts)/D))) *D + ts;
end
sync123
