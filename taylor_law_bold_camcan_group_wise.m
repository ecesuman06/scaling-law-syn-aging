clear all;clc;
close all

% load('camcan_fmri_data_analysis/rest_sorted.mat');tm=rest_sorted;strr1='resto5_you';strr2='resto5_mid';strr3='resto5_old';

% load('camcan_fmri_data_analysis/movie_sorted.mat');tm=movie_sorted;strr1='movieo5_you';strr2='movieo5_mid';strr3='movieo5_old';

load('camcan_fmri_data_analysis/smt_sorted.mat');tm=smt_sorted;strr1='smto5_you';strr2='smto5_mid';strr3='smto5_old';

age1=cell2mat(tm(:,3));
bold=tm(:,2);

sub_id=find(age1<=38); %%%%you_id;
sub = length(sub_id) ;
[b,p,R,sync,b_ana] = TL(bold(sub_id)) ;
noise_in = 1 ;
fsz = 28 ;
lfsz = 28 ;
% % % --------------------------------------------------------------------
[b_marg,p_marg,R_marg,sync_marg,b_ana_marg]=TL_surr(bold(sub_id));
b=b_ana;b_marg=b_ana_marg;
b_sync=b-b_marg;
% % % ------------------- sync vs b_sync -------------------------
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
    fill(points(1,:),points(2,:),'k','Facealpha',0.1)
    hold on
    plot(points(1,:),points(2,:),'r','linew',2)
    ylims = ylim; 
    ylim([0,ylims(2)]);
    axis tight
    axis off
    set(gca,'units','inches','position',[1.5 5 4 0.5])
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
    fill(points(2,:),points(1,:),'k','Facealpha',0.1)
    hold on
    plot(points(2,:),points(1,:),'r','linew',2)
%     ylims = ylim; 
%     ylim([0,ylims(2)]);
    set(gca,'units','inches','position',[5.5 1 0.5 4])
    axis tight
    axis off
subplot 223
    mdl1=fitlm(sync,data3);  %     anova(mdl)
    ab2=table2array(mdl1.Coefficients);
    s_vs_b=ab2(2,1);
    p_val=ab2(2,4);
    [rho,pval]=corr(sync,data3,'type','spearman');
    r2=mdl1.Rsquared.Ordinary;
    c = linspace(18,88,sub);
    scatter(sync,data3,50,[0.5 0.5 0.5],'filled')
%     cb=colorbar; cb.Location='south';  %'north' | 'south' | 'east' | 'west' | 'northoutside' | 'southoutside' | 'eastoutside' |'westoutside' | 'manual' | 'layout'
        hold on
    l=plot(mdl1);
    l(1).MarkerSize=0.1;
    l(2).LineWidth=2.3;
    l(2).Color='r';
%     l(2).Color=[135,62,35]/255;
    set(gca,'fontsize',fsz )
    legend off;    title ' ';
    set(gca,'xlim',[0,0.8],'ylim',[-1.25,0.2]);
    xlabel('sync, $\Omega$','interpreter','latex','Fontsize',lfsz)
    ylabel('slope, $b_{syn}$','interpreter','latex','Fontsize',lfsz)
%     str=sprintf('%s',num2str(p_val));
    str=sprintf('%s',num2str(round(pval,4)));    
    text(0.55,0.2,['$p$=',str],'Fontsize',fsz,'interpreter','latex')
    str=sprintf('%s',num2str(round(rho,3)));    
    text(0.4,0.2,['$\rho$=',str],'Fontsize',fsz,'interpreter','latex')
    str=sprintf('=%s',num2str(round(s_vs_b,4)));
    text(0.55,0.15,['$\beta$',str],'Fontsize',fsz,'interpreter','latex')
    str=sprintf('=%s',num2str(round(r2,4)));
    text(0.55,0.1,['$R^{2}$',str],'Fontsize',fsz,'interpreter','latex')
    set(gca,'units','inches','position',[1.5 1 4 4])
%     text(0.4,0.82,'age')
    axis tight
    box on
% saveas(gcf,['/home/mimi/Desktop/gaba_glut_new_analysis/camcan_fmri_data_analysis/figures/sync_b_sync_',strr1,'.fig']);
% saveas(gcf,['/home/mimi/Desktop/gaba_glut_new_analysis/camcan_fmri_data_analysis/figures/sync_b_sync_',strr1,'.png']);

% % % ==============================MIDDLE===============================
sub_id=find(45<=age1&age1<=60); %%%mid_id;
sub=length(sub_id);
[b,p,RR,sync,b_ana]=TL(bold(sub_id));
% % % --------------------------------------------------------------------
[b_marg,p_marg,RR_marg,sync_marg,b_ana_marg]=TL_surr(bold(sub_id));
b=b_ana;b_marg=b_ana_marg;
b_sync=b-b_marg;
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
    fill(points(1,:),points(2,:),'k','Facealpha',0.1)
    hold on
    plot(points(1,:),points(2,:),'r','linew',2)
    ylims = ylim; 
    ylim([0,ylims(2)]);
    axis tight
    axis off
    set(gca,'units','inches','position',[1.5 5 4 0.5])
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
    fill(points(2,:),points(1,:),'k','Facealpha',0.1)
    hold on
    plot(points(2,:),points(1,:),'r','linew',2)
%     ylims = ylim; 
%     ylim([0,ylims(2)]);
    set(gca,'units','inches','position',[5.5 1 0.5 4])
    axis tight
    axis off
subplot 223
    mdl1=fitlm(sync,data3);  %     anova(mdl)
    ab2=table2array(mdl1.Coefficients);
    s_vs_b=ab2(2,1);
    p_val=ab2(2,4);
    [rho,pval]=corr(sync,data3,'type','spearman');
    r2=mdl1.Rsquared.Ordinary;
    c = linspace(18,88,sub);
    scatter(sync,data3,50,[0.5 0.5 0.5],'filled')
%     cb=colorbar; cb.Location='south';  %'north' | 'south' | 'east' | 'west' | 'northoutside' | 'southoutside' | 'eastoutside' |'westoutside' | 'manual' | 'layout'
        hold on
    l=plot(mdl1);
    l(1).MarkerSize=0.1;
    l(2).LineWidth=2.3;
    l(2).Color='r';
%     l(2).Color=[135,62,35]/255;
    set(gca,'fontsize',fsz )
    legend off;    title ' ';
    set(gca,'xlim',[0,0.8],'ylim',[-1.25,0.2]);
    xlabel('sync, $\Omega$','interpreter','latex','Fontsize',lfsz)
    ylabel('slope, $b_{syn}$','interpreter','latex','Fontsize',lfsz)
%     str=sprintf('%s',num2str(p_val));
    str=sprintf('%s',num2str(round(pval,4)));    
    text(0.55,0.2,['$p$=',str],'Fontsize',fsz,'interpreter','latex')
    str=sprintf('%s',num2str(round(rho,3)));    
    text(0.4,0.2,['$\rho$=',str],'Fontsize',fsz,'interpreter','latex')
    str=sprintf('=%s',num2str(round(s_vs_b,4)));
    text(0.55,0.15,['$\beta$',str],'Fontsize',fsz,'interpreter','latex')
    str=sprintf('=%s',num2str(round(r2,4)));
    text(0.55,0.1,['$R^{2}$',str],'Fontsize',fsz,'interpreter','latex')
    set(gca,'units','inches','position',[1.5 1 4 4])
%     text(0.4,0.82,'age')
    axis tight
    box on
% saveas(gcf,['/home/mimi/Desktop/gaba_glut_new_analysis/camcan_fmri_data_analysis/figures/sync_b_sync_',strr2,'.fig']);
% saveas(gcf,['/home/mimi/Desktop/gaba_glut_new_analysis/camcan_fmri_data_analysis/figures/sync_b_sync_',strr2,'.png']);

% % % =========================OLD=================================
sub_id=find(65<=age1);%%%%%%%old
sub=length(sub_id);
[b,p,RR,sync,b_ana]=TL(bold(sub_id));
% % % --------------------------------------------------------------------
[b_marg,p_marg,RR_marg,sync_marg,b_ana_marg]=TL_surr(bold(sub_id));
b=b_ana;b_marg=b_ana_marg;
b_sync=b-b_marg;
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
    fill(points(1,:),points(2,:),'k','Facealpha',0.1)
    hold on
    plot(points(1,:),points(2,:),'r','linew',2)
    ylims = ylim; 
    ylim([0,ylims(2)]);
    axis tight
    axis off
    set(gca,'units','inches','position',[1.5 5 4 0.5])
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
    fill(points(2,:),points(1,:),'k','Facealpha',0.1)
    hold on
    plot(points(2,:),points(1,:),'r','linew',2)
%     ylims = ylim; 
%     ylim([0,ylims(2)]);
    set(gca,'units','inches','position',[5.5 1 0.5 4])
    axis tight
    axis off
subplot 223
    mdl1=fitlm(sync,data3);  %     anova(mdl)
    ab2=table2array(mdl1.Coefficients);
    s_vs_b=ab2(2,1);
    p_val=ab2(2,4);
    [rho,pval]=corr(sync,data3,'type','spearman');
    r2=mdl1.Rsquared.Ordinary;
    c = linspace(18,88,sub);
    scatter(sync,data3,50,[0.5 0.5 0.5],'filled')
%     cb=colorbar; cb.Location='south';  %'north' | 'south' | 'east' | 'west' | 'northoutside' | 'southoutside' | 'eastoutside' |'westoutside' | 'manual' | 'layout'
        hold on
    l=plot(mdl1);
    l(1).MarkerSize=0.1;
    l(2).LineWidth=2.3;
    l(2).Color='r';
%     l(2).Color=[135,62,35]/255;
    set(gca,'fontsize',fsz)
    legend off;    title ' ';
    set(gca,'xlim',[0,0.8],'ylim',[-1.25,0.2]);
    xlabel('sync, $\Omega$','interpreter','latex','Fontsize',lfsz)
    ylabel('slope, $b_{syn}$','interpreter','latex','Fontsize',lfsz)
%     str=sprintf('%s',num2str(p_val));
    str=sprintf('%s',num2str(round(pval,4)));    
    text(0.55,0.2,['$p$=',str],'Fontsize',fsz,'interpreter','latex')
    str=sprintf('%s',num2str(round(rho,3)));    
    text(0.4,0.2,['$\rho$=',str],'Fontsize',fsz,'interpreter','latex')
    str=sprintf('=%s',num2str(round(s_vs_b,4)));
    text(0.55,0.15,['$\beta$',str],'Fontsize',fsz,'interpreter','latex')
    str=sprintf('=%s',num2str(round(r2,4)));
    text(0.55,0.1,['$R^{2}$',str],'Fontsize',fsz,'interpreter','latex')
    set(gca,'units','inches','position',[1.5 1 4 4])
%     text(0.4,0.82,'age')
    axis tight
    box on

% saveas(gcf,['/home/mimi/Desktop/gaba_glut_new_analysis/camcan_fmri_data_analysis/figures/sync_b_sync_',strr3,'.fig']);
% saveas(gcf,['/home/mimi/Desktop/gaba_glut_new_analysis/camcan_fmri_data_analysis/figures/sync_b_sync_',strr3,'.png']);
