
p=[1,2,4,8,16];
t=[450.82,62.67,19.89,9.1,5.7];

speed_up =t(1)./t;

figure('color','w');

plot(p,t);


set(gca, 'YScale', 'log') ;

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


figure('color','w');

plot(p,speed_up);


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