clear;clc;
X=[1/(5-1), 1/(13-1), 1/(21-1), 1/(29-1), 1/(41-1)];

% three point to calculate RHS
Z=[1.0726033382850e-01, 1.2989016208354e-02,4.7101081268095e-03,2.4079632435519e-03,1.1811692077520e-03];

ld = 3; % linewidth
ft = 30;
figure('color','w');
plot(X,Z,'LineWidth',ld,'color','k');hold on;
scatter(X,Z,60,'filled');
set(gca, 'YScale', 'log') ;
set(gca, 'XScale', 'log') ;
% daspect([1 1 1]);

xlabel('h','FontSize',ft,'FontName','Times')
ylabel('||u-u_{h}||_{L2}','FontSize',ft,'FontName','Times')
set(gca,'FontSize',ft,'FontName','Times New Roman');


% plot triangle showing slope dx^2 for second order accuracy
xs1 = 0.07;
xs2 = 0.11;
ys1 =0.004;
ys2 = (xs2/xs1)^2 * ys1;
xs = [xs1,xs2,xs2,xs1];
ys = [ys1,ys1,ys2,ys1];
loglog(xs,ys,'r')

grid on
grid minor

axes = gca;

axes.MinorGridLineStyle = '-';
axes.MinorGridColor = [0.5 0.5 0.5];
axes.GridColor = [0.5 0.5 0.5];
axes.LineWidth = 1.1;


