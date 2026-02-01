clear all;clc;
close all
savefig=0;
load('rest_sorted.mat'); tm=rest_sorted; strr1='rest';
% load('movie_sorted.mat');tm=movie_sorted;strr1='movie';
% load( 'smt_sorted.mat' ) ; tm    = smt_sorted ; strr1 = 'smt' ;
age   = cell2mat(tm(:,3));
bold  = tm(:,2);  sub = size(bold,1);
[b,p,RR,sync,b_ana]   = TL(bold) ;
[b_marg,p_marg,RR_marg,sync_marg,b_marg_ana]=TL_surr(bold);
b_sync=b-b_marg;

% % % ------------------------Fitting----------------------------------
% % % % ---------------Sync vs b---------------------------------------
fontsiz=25;
latex_fontsiz=35;
figure
set(gcf,'units','inches','position',[2 2  6.2 5.8])

subplot 211
[heights,centers] = hist(sync);
n = length(centers);
w = centers(2)-centers(1);
t = linspace(centers(1)-w/2,centers(end)+w/2,n+1);
p = fix(n/2);
dt = diff(t);
Fvals = cumsum([0,heights.*dt]);
F = spline(t, [0, Fvals, 0]);
DF = fnder(F);  % computes its first derivative
% fnplt(DF, 'r')
[points] = fnplt(DF);
fill(points(1,:),points(2,:),'r','Facealpha',0.1)
hold on
plot(points(1,:),points(2,:),'r','linew',2)
ylims = ylim;
ylim([0,ylims(2)]);
%     axis tight
axis off
set(gca,'units','inches','position',[1.5 5.01 4 0.5],'xlim',[0,0.7])
%     title(strr1,'fontweight','normal')

subplot 224
data1=b;
[heights,centers] = hist(data1);
n = length(centers);
w = centers(2)-centers(1);
t = linspace(centers(1)-w/2,centers(end)+w/2,n+1);
p = fix(n/2);
dt = diff(t);
Fvals = cumsum([0,heights.*dt]);
F = spline(t, [0, Fvals, 0]);
DF = fnder(F);  % computes its first derivative
% fnplt(DF, 'r')
[points] = fnplt(DF);
fill(points(2,:),points(1,:),'r','Facealpha',0.1)
hold on
plot(points(2,:),points(1,:),'r','linew',2)
set(gca,'units','inches','position',[5.5 1 0.5 4],'ylim',[0,2.1])
%     axis tight
axis off

subplot 223
mdl1=fitlm(sync,data1);  %     anova(mdl)
ab2=table2array(mdl1.Coefficients);
s_vs_b=ab2(2,1);
p_val=ab2(2,4);
r2=mdl1.Rsquared.Ordinary;
c = linspace(18,88,sub);
scatter(sync,data1,50,c,'filled')
cb=colorbar;    cb.Location='south';  %'north' | 'south' | 'east' | 'west' | 'northoutside' | 'southoutside' | 'eastoutside' |'westoutside' | 'manual' | 'layout'
cb.Position=[0.27 0.24 0.4 0.016];
cb.FontSize= 15;
hold on
l=plot(mdl1);
l(1).MarkerSize=0.1;
l(2).LineWidth=2.3;
l(2).Color=[135,62,35]/255;
set(gca,'fontsize',fontsiz )
legend off
title ' '
xlabel('sync, $\Omega$','interpreter','latex','Fontsize',latex_fontsiz)
ylabel('slope, $b$','interpreter','latex','Fontsize',latex_fontsiz)
set(gca,'units','inches','position',[1.5 1 4 4],'xlim',[0,0.7],'ylim',[0,2.1])
text(0.16,0.35,'age','Fontsize',20)
box on
if savefig==1
saveas(gcf,['C:\Users\Mimi\Desktop\camcan_fmri_data_analysis\fig_age_b_age_sync\analytical/syn_b_',strr1,'_fit.png']);
end
% % % % % ---------------Sync vs b_marg---------------------------------------
figure
set(gcf,'units','inches','position',[2 2  6.2 5.8])
subplot 211
[heights,centers] = hist(sync);
n = length(centers);
w = centers(2)-centers(1);
t = linspace(centers(1)-w/2,centers(end)+w/2,n+1);
p = fix(n/2);
dt = diff(t);
Fvals = cumsum([0,heights.*dt]);
F = spline(t, [0, Fvals, 0]);
DF = fnder(F);  % computes its first derivative
% fnplt(DF, 'r')
[points] = fnplt(DF);
fill(points(1,:),points(2,:),'r','Facealpha',0.1)
hold on
plot(points(1,:),points(2,:),'r','linew',2)
ylims = ylim;
ylim([0,ylims(2)]);
axis off
set(gca,'units','inches','position',[1.5 5 4 0.5],'xlim',[0,0.7])
% title(strr1,'fontweight','normal')

subplot 224
data2=b_marg;
[heights,centers] = hist(data2);
n = length(centers);
w = centers(2)-centers(1);
t = linspace(centers(1)-w/2,centers(end)+w/2,n+1);
p = fix(n/2);
dt = diff(t);
Fvals = cumsum([0,heights.*dt]);
F = spline(t, [0, Fvals, 0]);
DF = fnder(F);  % computes its first derivative
% fnplt(DF, 'r')
[points] = fnplt(DF);
fill(points(2,:),points(1,:),'r','Facealpha',0.1)
hold on
plot(points(2,:),points(1,:),'r','linew',2)
set(gca,'units','inches','position',[5.5 1 0.5 4],'ylim',[0,2.1])
% axis tight
axis off

subplot 223
mdl1=fitlm(sync,data2);  %     anova(mdl)
ab2=table2array(mdl1.Coefficients);
s_vs_b=ab2(2,1);
p_val=ab2(2,4);
r2=mdl1.Rsquared.Ordinary;
c = linspace(18,88,sub);
scatter(sync,data2,50,c,'filled')
%     cb=colorbar; cb.Location='south';  %'north' | 'south' | 'east' | 'west' | 'northoutside' | 'southoutside' | 'eastoutside' |'westoutside' | 'manual' | 'layout'
hold on
l=plot(mdl1);
l(1).MarkerSize=0.1;
l(2).LineWidth=2.3;
l(2).Color=[135,62,35]/255;
set(gca,'fontsize',fontsiz )
legend off
title ' '
xlabel('sync, $\Omega$','interpreter','latex','Fontsize',latex_fontsiz)
ylabel('slope, $b_{marg}$','interpreter','latex','Fontsize',latex_fontsiz)
set(gca,'units','inches','position',[1.5 1 4 4],'xlim',[0,0.7],'ylim',[0,2.1])
%     axis tight;
box on
if savefig==1
saveas(gcf,['C:\Users\Mimi\Desktop\camcan_fmri_data_analysis\fig_age_b_age_sync\analytical/sync_bmarg_',strr1,'_fit.png']);
end
% % % ------------------- sync vs b_sync   -------------------------
figure
    set(gcf,'units','inches','position',[2 2  6.2 5.8])
subplot 211
    [heights,centers] = hist(sync);
    n = length(centers);
    w = centers(2)-centers(1);
    t = linspace(centers(1)-w/2,centers(end)+w/2,n+1);
    p = fix(n/2);
    dt = diff(t);
    Fvals = cumsum([0,heights.*dt]);
    F = spline(t, [0, Fvals, 0]);
    DF = fnder(F);  % computes its first derivative
    % fnplt(DF, 'r')
    [points] = fnplt(DF);
    fill(points(1,:),points(2,:),'r','Facealpha',0.1)
    hold on
    plot(points(1,:),points(2,:),'r','linew',2)
    ylims = ylim;
    ylim([0,ylims(2)]);
%     axis tight
    axis off
    set(gca,'units','inches','position',[1.5 5 4 0.5],'xlim',[0,0.7])
% title(strr1,'fontweight','normal')

subplot 224
    data3=b_sync;
    [heights,centers] = hist(data3);
    n = length(centers);
    w = centers(2)-centers(1);
    t = linspace(centers(1)-w/2,centers(end)+w/2,n+1);
    p = fix(n/2);
    dt = diff(t);
    Fvals = cumsum([0,heights.*dt]);
    F = spline(t, [0, Fvals, 0]);
    DF = fnder(F);  % computes its first derivative
    % fnplt(DF, 'r')
    [points] = fnplt(DF);
    fill(points(2,:),points(1,:),'r','Facealpha',0.1)
    hold on
    plot(points(2,:),points(1,:),'r','linew',2)
%     ylims = ylim;
%     ylim([0,ylims(2)]);
    set(gca,'units','inches','position',[5.5 1 0.5 4],'ylim',[-.7,0.7])
%     axis tight
    axis off

subplot 223
    mdl1=fitlm(sync,data3);  %     anova(mdl)
    ab2=table2array(mdl1.Coefficients);
    s_vs_b=ab2(2,1);
    p_val=ab2(2,4);
    [rho,pval]=corr(sync,data3,'type','spearman');
    r2=mdl1.Rsquared.Ordinary;
%     c = linspace(18,88,sub);
    scatter(sync,data3,50,age,'filled')
%     cb=colorbar; cb.Location='south';  %'north' | 'south' | 'east' | 'west' | 'northoutside' | 'southoutside' | 'eastoutside' |'westoutside' | 'manual' | 'layout'
    hold on
    l=plot(mdl1);
    l(1).MarkerSize=0.1;
    l(2).LineWidth=2.3;
    l(2).Color=[135,62,35]/255;
    set(gca,'fontsize',fontsiz )
    legend off;    title ' ';
    set(gca,'xlim',[0,0.8],'ylim',[-1.25,0.2]);
    xlabel('sync, $\Omega$','interpreter','latex','Fontsize',latex_fontsiz)
    ylabel('slope, $b_{syn}$','interpreter','latex','Fontsize',latex_fontsiz)
     set(gca,'units','inches','position',[1.5 1 4 4],'xlim',[0,0.7],'ylim',[-0.7,0.7])
%     str=sprintf('%s',num2str(p_val));
    str=sprintf('%s',num2str(round(p_val,3)));
    text(max(xlim)*0.65,max(ylim)*0.9,['$p$=',str],'Fontsize',fontsiz,'interpreter','latex')
    str=sprintf('%s',num2str(round(rho,3)));
    text(max(xlim)*0.65,max(ylim)*0.76,['$\rho$=',str],'Fontsize',fontsiz,'interpreter','latex')
    str=sprintf('=%s',num2str(round(s_vs_b,4)));
    text(max(xlim)*0.65,min(ylim)*0.76,['$\beta$',str],'Fontsize',fontsiz,'interpreter','latex')
    str=sprintf('=%s',num2str(round(r2,4)));
    text(max(xlim)*0.6,min(ylim)*0.9,['$R^{2}$',str],'Fontsize',fontsiz,'interpreter','latex')
   
%     text(0.4,0.82,'age')
%     axis tight
    box on
    if savefig==1
saveas(gcf,['C:\Users\Mimi\Desktop\camcan_fmri_data_analysis\fig_age_b_age_sync\analytical/sync_b_sync_',strr1,'_fit.png']);

    end
% % % %=========================ANALYTICAL============================% % % ==================================================================
b_sync_ana=b_ana-b_marg_ana;
% % % % ---------------Sync vs b---------------------------------------
fontsiz=22;
latex_fontsiz=35;
figure
set(gcf,'units','inches','position',[2 2  6.2 5.8])

subplot 211
[heights,centers] = hist(sync);
n = length(centers);
w = centers(2)-centers(1);
t = linspace(centers(1)-w/2,centers(end)+w/2,n+1);
p = fix(n/2);
dt = diff(t);
Fvals = cumsum([0,heights.*dt]);
F = spline(t, [0, Fvals, 0]);
DF = fnder(F);  % computes its first derivative
% fnplt(DF, 'r')
[points] = fnplt(DF);
fill(points(1,:),points(2,:),'r','Facealpha',0.1)
hold on
plot(points(1,:),points(2,:),'r','linew',2)
ylims = ylim;
ylim([0,ylims(2)]);
%     axis tight
axis off
set(gca,'units','inches','position',[1.5 5.01 4 0.5],'xlim',[0,0.7])
%     title(strr1,'fontweight','normal')

subplot 224
data1=b;
[heights,centers] = hist(data1);
n = length(centers);
w = centers(2)-centers(1);
t = linspace(centers(1)-w/2,centers(end)+w/2,n+1);
p = fix(n/2);
dt = diff(t);
Fvals = cumsum([0,heights.*dt]);
F = spline(t, [0, Fvals, 0]);
DF = fnder(F);  % computes its first derivative
% fnplt(DF, 'r')
[points] = fnplt(DF);
fill(points(2,:),points(1,:),'r','Facealpha',0.1)
hold on
plot(points(2,:),points(1,:),'r','linew',2)
set(gca,'units','inches','position',[5.5 1 0.5 4],'ylim',[0,2.1])
%     axis tight
axis off

subplot 223
mdl1=fitlm(sync,data1);  %     anova(mdl)
ab2=table2array(mdl1.Coefficients);
s_vs_b=ab2(2,1);
p_val=ab2(2,4);
r2=mdl1.Rsquared.Ordinary;
c = linspace(18,88,sub);
scatter(sync,data1,50,c,'filled')
cb=colorbar;    cb.Location='south';  %'north' | 'south' | 'east' | 'west' | 'northoutside' | 'southoutside' | 'eastoutside' |'westoutside' | 'manual' | 'layout'
cb.Position=[0.27 0.24 0.4 0.016];
cb.FontSize= 15;
hold on
l=plot(mdl1);
l(1).MarkerSize=0.1;
l(2).LineWidth=2.3;
l(2).Color=[135,62,35]/255;
set(gca,'fontsize',fontsiz )
legend off
title ' '
xlabel('sync, $\Omega$','interpreter','latex','Fontsize',latex_fontsiz)
ylabel('slope, $b$','interpreter','latex','Fontsize',latex_fontsiz)
set(gca,'units','inches','position',[1.5 1 4 4],'xlim',[0,0.7],'ylim',[0,2.1])
text(0.16,0.35,'age','Fontsize',20)
box on
if savefig==1
saveas(gcf,['C:\Users\Mimi\Desktop\camcan_fmri_data_analysis\fig_age_b_age_sync\analytical/syn_b_',strr1,'_ana.png']);
end
% % % % ---------------Sync vs b_marg---------------------------------------
figure
set(gcf,'units','inches','position',[2 2  6.2 5.8])
subplot 211
[heights,centers] = hist(sync);
n = length(centers);
w = centers(2)-centers(1);
t = linspace(centers(1)-w/2,centers(end)+w/2,n+1);
p = fix(n/2);
dt = diff(t);
Fvals = cumsum([0,heights.*dt]);
F = spline(t, [0, Fvals, 0]);
DF = fnder(F);  % computes its first derivative
% fnplt(DF, 'r')
[points] = fnplt(DF);
fill(points(1,:),points(2,:),'r','Facealpha',0.1)
hold on
plot(points(1,:),points(2,:),'r','linew',2)
ylims = ylim;
ylim([0,ylims(2)]);
axis off
set(gca,'units','inches','position',[1.5 5 4 0.5],'xlim',[0,0.7])
% title(strr1,'fontweight','normal')

subplot 224
data2=b_marg_ana;
[heights,centers] = hist(data2);
n = length(centers);
w = centers(2)-centers(1);
t = linspace(centers(1)-w/2,centers(end)+w/2,n+1);
p = fix(n/2);
dt = diff(t);
Fvals = cumsum([0,heights.*dt]);
F = spline(t, [0, Fvals, 0]);
DF = fnder(F);  % computes its first derivative
% fnplt(DF, 'r')
[points] = fnplt(DF);
fill(points(2,:),points(1,:),'r','Facealpha',0.1)
hold on
plot(points(2,:),points(1,:),'r','linew',2)
set(gca,'units','inches','position',[5.5 1 0.5 4],'ylim',[0,2.1])
% axis tight
axis off

subplot 223
mdl1=fitlm(sync,data2);  %     anova(mdl)
ab2=table2array(mdl1.Coefficients);
s_vs_b=ab2(2,1);
p_val=ab2(2,4);
r2=mdl1.Rsquared.Ordinary;
c = linspace(18,88,sub);
scatter(sync,data2,50,c,'filled')
%     cb=colorbar; cb.Location='south';  %'north' | 'south' | 'east' | 'west' | 'northoutside' | 'southoutside' | 'eastoutside' |'westoutside' | 'manual' | 'layout'
hold on
l=plot(mdl1);
l(1).MarkerSize=0.1;
l(2).LineWidth=2.3;
l(2).Color=[135,62,35]/255;
set(gca,'fontsize',fontsiz )
legend off
title ' '
xlabel('sync, $\Omega$','interpreter','latex','Fontsize',latex_fontsiz)
ylabel('slope, $b_{marg}$','interpreter','latex','Fontsize',latex_fontsiz)
set(gca,'units','inches','position',[1.5 1 4 4],'xlim',[0,0.7],'ylim',[0,2.1])
%     axis tight;
box on
% saveas(gcf,['/home/mimi/Desktop/gaba_glut_new_analysis/camcan_fmri_data_analysis/fig_age_b_age_sync/sync_bmarg_',strr1,'.fig']);
if savefig==1
saveas(gcf,['C:\Users\Mimi\Desktop\camcan_fmri_data_analysis\fig_age_b_age_sync\analytical/sync_bmarg_',strr1,'_ana.png']);
end
% % % %======================================================================
% % % ------------------- sync vs b_sync   -------------------------
figure
    set(gcf,'units','inches','position',[2 2  6.2 5.8])
subplot 211
    [heights,centers] = hist(sync);
    n = length(centers);
    w = centers(2)-centers(1);
    t = linspace(centers(1)-w/2,centers(end)+w/2,n+1);
    p = fix(n/2);
    dt = diff(t);
    Fvals = cumsum([0,heights.*dt]);
    F = spline(t, [0, Fvals, 0]);
    DF = fnder(F);  % computes its first derivative
    % fnplt(DF, 'r')
    [points] = fnplt(DF);
    fill(points(1,:),points(2,:),'r','Facealpha',0.1)
    hold on
    plot(points(1,:),points(2,:),'r','linew',2)
    ylims = ylim;
    ylim([0,ylims(2)]);
%     axis tight
    axis off
    set(gca,'units','inches','position',[1.5 5 4 0.5],'xlim',[0,0.7])
% title(strr1,'fontweight','normal')

subplot 224
    data3=b_sync_ana;
    [heights,centers] = hist(data3);
    n = length(centers);
    w = centers(2)-centers(1);
    t = linspace(centers(1)-w/2,centers(end)+w/2,n+1);
    p = fix(n/2);
    dt = diff(t);
    Fvals = cumsum([0,heights.*dt]);
    F = spline(t, [0, Fvals, 0]);
    DF = fnder(F);  % computes its first derivative
    % fnplt(DF, 'r')
    [points] = fnplt(DF);
    fill(points(2,:),points(1,:),'r','Facealpha',0.1)
    hold on
    plot(points(2,:),points(1,:),'r','linew',2)
%     ylims = ylim;
%     ylim([0,ylims(2)]);
    set(gca,'units','inches','position',[5.5 1 0.5 4],'ylim',[-.7,0.7])
%     axis tight
    axis off

subplot 223
    mdl1=fitlm(sync,data3);  %     anova(mdl)
    ab2=table2array(mdl1.Coefficients);
    s_vs_b=ab2(2,1);
    p_val=ab2(2,4);
    [rho,pval]=corr(sync,data3,'type','spearman');
    r2=mdl1.Rsquared.Ordinary;
%     c = linspace(18,88,sub);
    scatter(sync,data3,50,age,'filled')
%     cb=colorbar; cb.Location='south';  %'north' | 'south' | 'east' | 'west' | 'northoutside' | 'southoutside' | 'eastoutside' |'westoutside' | 'manual' | 'layout'
    hold on
    l=plot(mdl1);
    l(1).MarkerSize=0.1;
    l(2).LineWidth=2.3;
    l(2).Color=[135,62,35]/255;
    set(gca,'fontsize',fontsiz )
    legend off;    title ' ';
    xlabel('sync, $\Omega$','interpreter','latex','Fontsize',latex_fontsiz)
    ylabel('slope, $b_{syn}$','interpreter','latex','Fontsize',latex_fontsiz)
     set(gca,'units','inches','position',[1.5 1 4 4],'xlim',[0,0.7],'ylim',[-0.7,0.7])
%     str=sprintf('%s',num2str(p_val));
    str=sprintf('%s',num2str(round(p_val,3)));
    text(max(xlim)*0.65,max(ylim)*0.9,['$p$=',str],'Fontsize',fontsiz,'interpreter','latex')
    str=sprintf('%s',num2str(round(rho,3)));
    text(max(xlim)*0.65,max(ylim)*0.76,['$\rho$=',str],'Fontsize',fontsiz,'interpreter','latex')
    str=sprintf('=%s',num2str(round(s_vs_b,4)));
    text(max(xlim)*0.65,min(ylim)*0.76,['$\beta$',str],'Fontsize',fontsiz,'interpreter','latex')
    str=sprintf('=%s',num2str(round(r2,4)));
    text(max(xlim)*0.6,min(ylim)*0.9,['$R^{2}$',str],'Fontsize',fontsiz,'interpreter','latex')
   
%     text(0.4,0.82,'age')
%     axis tight
    box on
% saveas(gcf,['/home/mimi/Desktop/gaba_glut_new_analysis/camcan_fmri_data_analysis/fig_age_b_age_sync/sync_b_sync_',strr1,'.fig']);
if savefig==1
saveas(gcf,['C:\Users\Mimi\Desktop\camcan_fmri_data_analysis\fig_age_b_age_sync\analytical/sync_b_sync_',strr1,'_ana.png']);
end












% % % ===================================================================
% figure;
% set(gcf,'units','inches','position',[2 2  18 6])
% tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact')
% % % % =============== age_vs_b ========================================
% nexttile
%   ylabels=['slope, $b$'];
%   plot_age_b_bsyn_syn(age,b,ylabels)
% % ==================== age_vs_bsync ================================
% nexttile
%   ylabels=['slope, $b_{syn}$'];
%     plot_age_b_bsyn_syn(age,b_sync,ylabels)
% % % % =============== age_vs_sync ========================================
% nexttile
%     ylabels=['sync, $\Omega$'];
%     plot_age_b_bsyn_syn(age,sync,ylabels)
% %     saveas(gcf,['/home/mimi/Desktop/gaba_glut_new_analysis/camcan_fmri_data_analysis/fig_age_b_age_sync/age_b_bsyn_syn_',strr1,'.fig'])
% %     saveas(gcf,['/home/mimi/Desktop/gaba_glut_new_analysis/camcan_fmri_data_analysis/fig_age_b_age_sync/age_b_bsyn_syn_',strr1,'.png'])
%
%
