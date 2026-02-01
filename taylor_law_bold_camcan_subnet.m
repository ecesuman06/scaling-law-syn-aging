clear all;clc;
close all

% load('rest_sorted.mat');tm=rest_sorted;strr1='resto5';
load('movie_sorted.mat'); tm=movie_sorted; strr1='movieo5';
% load('smt_sorted.mat'); tm=smt_sorted; strr1='smto5';
age=cell2mat(tm(:,3));
bold=tm(:,2);sub=size(bold,1); 

% bm=cell2mat(bold(641));
% bm(:,109)=[];
% figure
% plot(bm)
% axis tight
% xlabel('time')
% ylabel('amplitude')
% set(gca,'fontsize',18)
% 
% figure
% histogram(bm(:,:),100)
% xlabel('values')
% ylabel('counts')
% set(gca,'fontsize',18)
% 
% bm=detrend(bm);
% figure
% plot(bm)
% axis tight
% xlabel('time')
% ylabel('amplitude')
% set(gca,'fontsize',18)
% 
% figure
% histogram(bm(:,:),100)
% xlabel('values')
% ylabel('counts')
% set(gca,'fontsize',18)


% subnet = [9, 10, 15,16,31,32,35,36,65,66,67,68]; str2='DMN';
% subnet=[7,8,13,14,59,60,61,62,65,66]; str2='FPN';
% subnet=[1,2,17,18,19,20,57,58,81,82]; str2='SMN';
subnet=[43:56];   str2='VN';
% subnet=[31:40,83,84,87,88];  str2='limbic';
% subnet=[41,42,71:78];  str2='subcort';
% subnet=[91:108]; str2='cerebellum';

[b,p,R,sync,b_ana]=TL_subnet(bold,subnet);
% % % --------------------------------------------------------------------
[b_marg,p_marg,R_marg,sync_marg,b_ana_marg]=TL_surr_subnet(bold,subnet);
b=b_ana;b_marg=b_ana_marg;
b_sync=b-b_marg;
% % % --------------------------------------------------------------------
% % % % % ---------------Sync vs b---------------------------------------
fontsiz=26;
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
    axis tight
    axis off
    set(gca,'units','inches','position',[1.5 5 4 0.5])
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
    set(gca,'units','inches','position',[5.5 1 0.5 4])
    axis tight
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
    cb.FontSize= 12;
        hold on
    l=plot(mdl1);
    l(1).MarkerSize=0.1;
    l(2).LineWidth=2.3;
    l(2).Color=[135,62,35]/255;
    set(gca,'fontsize',fontsiz )
    legend off
    title ' '
    set(gca,'xlim',[0,0.8])
    set(gca,'ylim',[-1.25,0.2])
    xlabel('sync, $\Omega$','interpreter','latex','Fontsize',latex_fontsiz)
    ylabel('slope, $b$','interpreter','latex','Fontsize',latex_fontsiz)
%     str=sprintf('%s',num2str(p_val));
%     text(0.48,2,['$p$=',str],'Fontsize',16,'interpreter','latex')
    str=sprintf('=%s',num2str(round(s_vs_b,4)));
    text(0.48,1.8,['$\beta$',str],'Fontsize',fontsiz,'interpreter','latex')
    str=sprintf('=%s',num2str(round(r2,4)));
    text(0.48,1.7,['$R^{2}$',str],'Fontsize',fontsiz,'interpreter','latex')
    set(gca,'units','inches','position',[1.5 1 4 4])
    text(0.275,0.45,'age','Fontsize',fontsiz)
    axis tight;box on
% saveas(gcf,['/home/mimi/Desktop/gaba_glut_new_analysis/camcan_fmri_data_analysis/figures/syn_b_',strr1,'.fig']);
% saveas(gcf,['/home/mimi/Desktop/gaba_glut_new_analysis/camcan_fmri_data_analysis/figures/syn_b_',strr1,'.png']);

% % %======================================================================
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
    axis tight
    axis off
    set(gca,'units','inches','position',[1.5 5 4 0.5])
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
    set(gca,'units','inches','position',[5.5 1 0.5 4])
    axis tight
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
    set(gca,'xlim',[0,0.8])
    set(gca,'ylim',[-1.25,0.2])
    xlabel('sync, $\Omega$','interpreter','latex','Fontsize',latex_fontsiz)
    ylabel('slope, $b_{marg}$','interpreter','latex','Fontsize',latex_fontsiz)
%     str=sprintf('%s',num2str(p_val));
%     text(0.48,2,['$p$=',str],'Fontsize',16,'interpreter','latex')
    str=sprintf('=%s',num2str(round(s_vs_b,4)));
    text(0.48,1.8,['$\beta$',str],'Fontsize',fontsiz,'interpreter','latex')
    str=sprintf('=%s',num2str(round(r2,4)));
    text(0.48,1.7,['$R^{2}$',str],'Fontsize',fontsiz,'interpreter','latex')
    set(gca,'units','inches','position',[1.5 1 4 4])
%     text(0.4,0.82,'age','Fontsize',14)
    axis tight;box on
% saveas(gcf,['/home/mimi/Desktop/gaba_glut_new_analysis/camcan_fmri_data_analysis/figures/sync_bmarg_',strr1,'.fig']);
% saveas(gcf,['/home/mimi/Desktop/gaba_glut_new_analysis/camcan_fmri_data_analysis/figures/sync_bmarg_',strr1,'.png']);

% % %======================================================================
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
    axis tight
    axis off
    set(gca,'units','inches','position',[1.5 5.2 4 0.5])
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
    set(gca,'units','inches','position',[5.5 1.2 0.5 4])
    axis tight
    axis off
subplot 223
   set(gca,'units','inches','position',[1.5 1.2 4 4])
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
%     str=sprintf('%s',num2str(p_val));
    str=sprintf('%s',num2str(round(pval,4)));    
    text(0.45,0.15,['$p$=',str],'Fontsize',fontsiz-1,'interpreter','latex')
    str=sprintf('%s',num2str(round(rho,3)));    
    text(0.45,0.12,['$\rho$=',str],'Fontsize',fontsiz-1,'interpreter','latex')
    str=sprintf('=%s',num2str(round(s_vs_b,4)));
    text(0.45,0.09,['$\beta$',str],'Fontsize',fontsiz-1,'interpreter','latex')
    str=sprintf('=%s',num2str(round(r2,4)));
    text(0.45,0.0,['$R^{2}$',str],'Fontsize',fontsiz-1,'interpreter','latex')
 %     text(0.4,0.82,'age')
    axis tight
    box on
ax=gca;
ax.LineWidth=2;

% saveas(gcf,['figures/sync_b_sync_',strr1,'-',str2,'.fig']);

