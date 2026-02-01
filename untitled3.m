clear all;  clc;  close all

load('camcan_fmri_data_analysis/rest_sorted.mat');
tm=rest_sorted;
strr1='rest';
bold=tm(:,2);sub=size(bold,1); 
bm=cell2mat(bold(641));
bm(:,109)=[];

figure
plot(bm)
% figure
% histogram(bm(:,:),100)


s1=phase_rand_surr(bm,1);

figure
plot(s1)
% figure
% histogram(s1(:,:),100)


s2 = PhaseRand_surrogates(bm, 1);
figure
plot(s2)
% figure
% histogram(s2(:,:),100)

figure
plot(bm(:,101),'k')
hold
plot(s1(:,101),'r')
plot(s2(:,101),'b')