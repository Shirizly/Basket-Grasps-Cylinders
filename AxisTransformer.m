% open figure of U
fig = gcf;
axObjs = fig.Children;
axObjs2 = axObjs(2);
dataObjs2 = axObjs2.Children;
x = dataObjs2(1).XData
y = dataObjs2(1).YData
z = dataObjs2(1).ZData
%% varying theta 1
subFolderName = '\Depth_r_05_h_2';
y_new = 2*(pi()-y);
zt = z.';
fig5 = figure;
[~,cf] = contourf(x,y_new,z,30);
% set(cf,'LineColor','none')
cb = colorbar;
xlabel('$\theta_1$','interpreter','latex','FontSize',20)
ylabel('arc-length','interpreter','latex','FontSize',20)
set(get(cb,'label'),'string','$\:U$','interpreter','latex','FontSize',20,'rotation',0);
set(gca,'YTick',[0.5*pi 0.6*pi 0.7*pi 0.8*pi 0.9*pi],'YTickLabel',{'0.5$\pi r$','0.6$\pi r$','0.7$\pi r$','0.8$\pi r$','0.9$\pi r$'},'tickLabelInterpreter','latex')
xlim([min(x),max(x(zt(:,1)>0))])
savefig(fig5,[pwd '\' subFolderName '\U_NewPhi'])
saveas(fig5,[pwd '\' subFolderName '\U_NewPhi.eps'],'epsc')
saveas(fig5,[pwd '\' subFolderName '\U_NewPhi.png'],'png')
%% varying ds
subFolderName = '\Depth_r_05_h_16_ds';
y_new = 2*(pi()-y);
x_new = x/0.8;
zt = z.';
fig5 = figure;
[~,cf] = contourf(x_new,y_new,z,30);
% set(cf,'LineColor','none')
cb = colorbar;
xlabel('$\delta s''$','interpreter','latex','FontSize',20)
ylabel('arc-length','interpreter','latex','FontSize',20)
set(get(cb,'label'),'string','$\:U$','interpreter','latex','FontSize',20,'rotation',0);
set(gca,'YTick',[0.3*pi 0.4*pi 0.5*pi 0.6*pi 0.7*pi 0.8*pi 0.9*pi],'YTickLabel',{'0.3$\pi r$','0.4$\pi r$','0.5$\pi r$','0.6$\pi r$','0.7$\pi r$','0.8$\pi r$','0.9$\pi r$'},'tickLabelInterpreter','latex')
xlim([min(x_new),max(x_new(zt(:,1)>0))])
ylims = get(gca,'ylim');
ylims(1) = 0.3*pi();
ylim(ylims);
savefig(fig5,[pwd '\' subFolderName '\U_NewPhi'])
saveas(fig5,[pwd '\' subFolderName '\U_NewPhi.eps'],'epsc')
saveas(fig5,[pwd '\' subFolderName '\U_NewPhi.png'],'png')
