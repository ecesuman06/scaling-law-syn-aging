clear all; clc; 
close all
% delete(gcp); c = parcluster;  parpool(c);
% rest_sorted = importdata('camcan_fmri_data_analysis/rest_sorted.mat');
% bold        = rest_sorted(:,2);
% age1    = cell2mat(rest_sorted(:,3));
% sub_id  = find(age1<=38);            %%%------you_id;
% sub_id2 = find(45<=age1&age1<=60);   %%%------mid_id;
% sub_id3 = find(65<=age1);            %%%------old_id;
% smt_sorted   = importdata('camcan_fmri_data_analysis/smt_sorted.mat');
% bold2        = smt_sorted(:,2);
% movie_sorted = importdata('camcan_fmri_data_analysis/movie_sorted.mat');
% bold3        = movie_sorted(:,2) ;
% noise_in = 1 ;
% parfor i = 1:1024   
%     disp( num2str(i) )
%     % % % ======================YOU=============================  
%     [ b, p, R, sync, b_ana ] = TL_2( bold( sub_id ) ) ;
%     [ b_marg, p_marg, R_marg, sync_marg, b_marg_ana ] = taylor_law_bold_camcan_surr_2( bold(sub_id), noise_in ) ;
%     b_sync     = b - b_marg ;
%     b_sync_ana = b_ana - b_marg_ana ;
%     %     figure
%     %     set(gcf,'units','inches','position',[2 2 15 18])
%     %     tiledlayout(3,3,'TileSpacing','Compact','Padding','Compact')
%     %     nexttile;
%     beta1_rest     = plot_syn_bsyn( sync, b_sync ) ;
%     beta1_rest_ana = plot_syn_bsyn( sync, b_sync_ana ) ;
%     % % % ======================MID=============================  
%     [ b2, p, RR, sync2, b2_ana ] = TL_2( bold(sub_id2) ) ;
%     [ b2_marg, p_marg, RR_marg, sync_marg, b2_marg_ana ] = taylor_law_bold_camcan_surr_2( bold(sub_id2), noise_in ) ;
%     b_sync2     = b2     - b2_marg ;
%     b_sync2_ana = b2_ana - b2_marg_ana;
%     %     nexttile;
%     beta2_rest      = plot_syn_bsyn( sync2, b_sync2 ) ;
%     beta2_rest_ana  = plot_syn_bsyn( sync2, b_sync2_ana ) ;
%     % % % =======================OLD=================================
%     [ b3, p, RR, sync3, b3_ana ] = TL_2(bold(sub_id3));
%     [ b3_marg, p_marg, RR_marg, sync_marg, b3_marg_ana ] = taylor_law_bold_camcan_surr_2(bold(sub_id3),noise_in);
%     b_sync3     = b3 - b3_marg ;
%     b_sync3_ana = b3_ana - b3_marg_ana ;
%     %     nexttile;
%     beta3_rest     = plot_syn_bsyn( sync3, b_sync3 ) ;
%     beta3_rest_ana = plot_syn_bsyn( sync3, b_sync3_ana ) ;
%     
%     % % % %---------------SMT--------------------------------
%     % % % ================ YOU ===================     
%     [b, p, R, sync, b_ana] = TL_2( bold2( sub_id ) ) ;
%     [b_marg, p_marg, R_marg, sync_marg, b_marg_ana] = taylor_law_bold_camcan_surr_2( bold2( sub_id ), noise_in ) ;
%     b_sync_smt      = b - b_marg  ;
%     b_sync_smt_ana  = b_ana-b_marg_ana ;
%     %     nexttile;
%     beta1_smt       = plot_syn_bsyn(sync,b_sync_smt);
%     beta1_smt_ana   = plot_syn_bsyn(sync,b_sync_smt_ana);
%     % % % ================= MID ======================
%     [b2, p, RR, sync2, b2_ana] = TL_2( bold2( sub_id2 ) ) ;
%     [b2_marg, p_marg, RR_marg, sync_marg, b2_marg_ana] = taylor_law_bold_camcan_surr_2( bold2(sub_id2), noise_in ) ;
%     b2_sync_smt     = b2 - b2_marg ;
%     b2_sync_smt_ana = b2_ana - b2_marg_ana ;
%     %     nexttile;
%     beta2_smt       = plot_syn_bsyn( sync2, b2_sync_smt);
%     beta2_smt_ana   = plot_syn_bsyn( sync2, b2_sync_smt_ana);
%     % % % ================ OLD ======================
%     [b3, p, RR, sync3, b3_ana ] = TL_2( bold2(sub_id3) ) ;
%     [b3_marg, p_marg, RR_marg, sync_marg, b3_marg_ana ] = taylor_law_bold_camcan_surr_2( bold2(sub_id3), noise_in ) ;
%     b3_sync_smt     = b3 - b3_marg ;
%     b3_sync_smt_ana = b3_ana - b3_marg_ana ;
%     %     nexttile;
%     beta3_smt       = plot_syn_bsyn( sync3, b3_sync_smt ) ;
%     beta3_smt_ana   = plot_syn_bsyn( sync3, b3_sync_smt_ana ) ;
%     
%     % % % %---------------MOVIE--------------------------------
%     %%% ================ YOU ===================    
%     [b, p, R, sync, b_ana] = TL_2( bold3( sub_id ) );
%     [b_marg, p_marg, R_marg, sync_marg, b_marg_ana] = taylor_law_bold_camcan_surr_2( bold3(sub_id), noise_in ) ;
%     b_sync_mov         = b - b_marg;
%     b_sync_mov_ana     = b_ana - b_marg_ana;
%     %     nexttile;
%     beta1_mov          = plot_syn_bsyn( sync, b_sync_mov ) ;
%     beta1_mov_ana      = plot_syn_bsyn( sync, b_sync_mov_ana ) ;
%     % % % ============== MID =====================
%     [ b2, p, RR, sync2, b2_ana ]  =  TL_2( bold3( sub_id2 ) ) ;
%     [ b2_marg, p_marg, RR_marg, sync_marg, b2_marg_ana ] = taylor_law_bold_camcan_surr_2( bold3(sub_id2), noise_in ) ;
%     b2_sync_mov       = b2 - b2_marg ;
%     b2_sync_mov_ana   = b2_ana - b2_marg_ana ;
%     %     nexttile;
%     beta2_mov         = plot_syn_bsyn( sync2, b2_sync_mov ) ;
%     beta2_mov_ana     = plot_syn_bsyn( sync2, b2_sync_mov_ana ) ;
%     % % % ===============OLD=====================
%     [ b3, p, RR, sync3, b3_ana ] = TL_2(bold3(sub_id3));
%     [b3_marg, p_marg, RR_marg, sync_marg, b3_marg_ana ] = taylor_law_bold_camcan_surr_2(bold3(sub_id3),noise_in);
%     b3_sync_mov       = b3     - b3_marg ;
%     b3_sync_mov_ana   = b3_ana - b3_marg_ana ;
%     %     nexttile;
%     beta3_mov         = plot_syn_bsyn( sync3, b3_sync_mov ) ;
%     beta3_mov_ana     = plot_syn_bsyn( sync3, b3_sync_mov_ana ) ;
%     
%     b1_fit = [ beta1_rest, beta2_rest, beta3_rest, beta1_smt, beta2_smt, beta3_smt, beta1_mov, beta2_mov, beta3_mov ] ;
%     b1_ana = [ beta1_rest_ana, beta2_rest_ana, beta3_rest_ana, beta1_smt_ana, beta2_smt_ana, beta3_smt_ana, beta1_mov_ana, beta2_mov_ana, beta3_mov_ana ] ;
%     
%     fi=fopen('C:\Users\Mimi\Desktop\camcan_fmri_data_analysis\trial_wise_b_rest_smt_mov_fit_1000run4.txt','a');
% %     "C:\Users\Mimi\Desktop\camcan_fmri_data_analysis\trial_wise_b_rest_smt_mov.txt"
%     fprintf(fi,'%s ',num2str(b1_fit));
%     fprintf(fi,'\n');
%     fclose(fi);
%     
%     fi1=fopen('C:\Users\Mimi\Desktop\camcan_fmri_data_analysis\trial_wise_b_rest_smt_mov_ana_1000run4.txt','a');
%     fprintf(fi1,'%s ',num2str(b1_ana));
%     fprintf(fi1,'\n');
%     fclose(fi1);
%     
% end

% % ________Result without standarizing BOLD ts dividing by std_________
% beta=importdata('camcan_fmri_data_analysis/trial_wise_b_rest_smt_mov.txt');
beta=importdata('C:\Users\Mimi\Desktop\camcan_fmri_data_analysis\trial_wise_b_rest_smt_mov_fit_1000run4.txt');
mbeta=mean(beta);
sbeta=std(beta);
figure
plot(1:3,mbeta(1:3),'-ro','markersize',6,'linew',2,'Markerfacecolor','r');;hold on
plot(1:3,mbeta(4:6),'-bo','markersize',6,'markerfacecolor','b','linew',2);;hold on
plot(1:3,mbeta(7:9),'-ko','markersize',6,'markerfacecolor','k','linew',2);;hold on
errorbar(1:3,mbeta(1:3),sbeta(1:3),'.r','linew',2);hold on
errorbar(1:3,mbeta(4:6),sbeta(4:6),'.b','linew',2)
errorbar(1:3,mbeta(7:9),sbeta(7:9),'.k','linew',2)
l=legend('rest','smt', 'movie');
l.EdgeColor='w';
l.FontSize=17;
set(gca,'xlim',[0.5,3.5],'ylim',[0,1],'fontsize',25)
set(gca,'xtick',1:3,'xticklabel',{'YA', 'MA', 'OA'})
ylabel('$\beta$','interpreter','latex','fontsize',35,'fontweight','bold')
%
%
% saveas(gcf,'C:\Users\Mimi\Desktop\camcan_fmri_data_analysis\fig_age_b_age_sync\analytical/group_b_rest_smt_mov_fit_1000run4.png')

beta=importdata('C:\Users\Mimi\Desktop\camcan_fmri_data_analysis\trial_wise_b_rest_smt_mov_ana_1000run4.txt');
mbeta=mean(beta);
sbeta=std(beta);
figure
plot(1:3,mbeta(1:3),'-ro','markersize',6,'linew',2,'Markerfacecolor','r'); hold on
plot(1:3,mbeta(4:6),'-bo','markersize',6,'markerfacecolor','b','linew',2); hold on
plot(1:3,mbeta(7:9),'-ko','markersize',6,'markerfacecolor','k','linew',2); hold on
errorbar(1:3,mbeta(1:3),sbeta(1:3),'.r','linew',2);hold on
errorbar(1:3,mbeta(4:6),sbeta(4:6),'.b','linew',2)
errorbar(1:3,mbeta(7:9),sbeta(7:9),'.k','linew',2)
% l=legend('rest','smt', 'movie');
% l.EdgeColor='w';
% l.FontSize=17;
set(gca,'xlim',[0.5,3.5],'ylim',[0,1],'fontsize',25)
set(gca,'xtick',1:3,'xticklabel',{'YA', 'MA', 'OA'})
ylabel('$\beta$','interpreter','latex','fontsize',35,'fontweight','bold')
% saveas(gcf,'C:\Users\Mimi\Desktop\camcan_fmri_data_analysis\fig_age_b_age_sync\fitted_results/group_b_rest_smt_mov_ana_1000run4.png')
