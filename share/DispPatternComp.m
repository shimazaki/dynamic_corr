function DispPattern(emp_p,model_p,c)

p = loglog(emp_p,model_p,'.');
set(p,'Color',c, 'MarkerSize', 3);

set(gca,'XScale','log','YScale','log');
xlabel('empirical'); ylabel('model');
axis square;
