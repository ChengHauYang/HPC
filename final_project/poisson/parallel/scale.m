
p=[1,2,4,8,16];
% t=[874.79,144.16,41.59,14.28,6.02]; % 161
t=[20.6,5.54,1.75,0.74,0.42]; 
Size =20;

speed_up =t(1)./t;

figure('color','w');

plot(p,t);hold on;
scatter(p,t,'filled');

set(gca, 'YScale', 'log') ;
set(gca, 'XScale', 'log') ;


set(gca,'FontSize',Size,'FontName','Times New Roman');
grid on
grid minor

axes = gca;

axes.MinorGridLineStyle = '-';
axes.MinorGridColor = [0.5 0.5 0.5];
axes.GridColor = [0.5 0.5 0.5];
axes.LineWidth = 1.1;


xlabel('Number of Processors');
ylabel('Time');

hold on;
% plot triangle showing slope dx^2 for second order accuracy
xs1 = 3;
xs2 = 5;
ys1 =6;
ys2 = ys1/ (xs2/xs1)^2;
xs = [xs1,xs2,xs2,xs1];
ys = [ys1,ys1,ys2,ys1];
loglog(xs,ys,'r')
daspect([1 1 1]);

%%

figure('color','w');

plot(p,speed_up);hold on;
scatter(p,speed_up,'filled');

% set(gca, 'YScale', 'log') ;

set(gca,'FontSize',Size,'FontName','Times New Roman');
grid on
grid minor

axes = gca;

axes.MinorGridLineStyle = '-';
axes.MinorGridColor = [0.5 0.5 0.5];
axes.GridColor = [0.5 0.5 0.5];
axes.LineWidth = 1.1;


xlabel('Number of Processors');
ylabel('Speedup');