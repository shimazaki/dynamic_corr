function DispHyperparam(hyper)

Q = hyper.Q;
F = hyper.F;
[M,M] = size(Q);

%%%%%% Covariance Matrix %%%%%%%%%
PCOLOR = zeros(M+1,M+1); 
subplot(2,1,1); PCOLOR(1:end-1,1:end-1) = Q-diag(diag(Q)); pcolor(PCOLOR); 
axis ij; axis square; colorbar; title('Q'); axis off; shading flat;
subplot(2,1,2); PCOLOR(1:end-1,1:end-1) = F-diag(diag(F)); pcolor(PCOLOR); 
axis ij; axis square; colorbar; title('F'); axis off; shading flat;

