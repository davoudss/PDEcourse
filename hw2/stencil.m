clc; clear all; clf

x = [1 0 -1 0 1 -1];
y = [0 1 0 -1 -1 1];
figure(1)
hold on
for k= 1:numel(x)
    plot([x(k) 0],[y(k) 0],'r-o','MarkerFaceColor','r','markersize',4)
    axis([-2 2 -2 2])
end
grid on


set(gca,'xtick',[-2 -1 0 1 2])
set(gca,'ytick',[-2 -1 0 1 2])
xtick = ['i-2';'i-1';' i ';'i+1';'i+2'];
ytick = ['j-2';'j-1';' j ';'j+1';'j+2'];
set(gca,'XTickLabel',xtick)
set(gca,'YTickLabel',ytick)
publication_fig
box on

figure(2)
x = [-1 1];
y = [-1 1];
hold on
for k= 1:numel(x)
    plot([x(k) 0],[y(k) 0],'r-o','MarkerFaceColor','r','markersize',4)
    axis([-2 2 -2 2])
end
grid on


set(gca,'xtick',[-2 -1 0 1 2])
set(gca,'ytick',[-2 -1 0 1 2])
xtick = ['i-2';'i-1';' i ';'i+1';'i+2'];
ytick = ['j-2';'j-1';' j ';'j+1';'j+2'];
set(gca,'XTickLabel',xtick)
set(gca,'YTickLabel',ytick)
publication_fig
box on