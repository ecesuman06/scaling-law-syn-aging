% function [b_fit,p_val,r2,sync_fc,b_analytical]=TL_nki(bold)
% % % % bold : input data of size: time x location
% sub=size(bold,1);  % total number of subjects
% p_val=zeros(sub,1); r2=zeros(sub,1);
% b_fit=zeros(sub,1);  b_actual=zeros(sub,1); b_analytical=zeros(sub,1);
% sync_rho=zeros(sub,1); sync_mu=zeros(sub,1); sync_fc=zeros(sub,1);
% % var_m=zeros(sub,1); var_mu=zeros(sub,1);
% 
% for i=1:sub
%     s1 = [] ;
%     s1 = bold{i} ;
% %     s1=s1./std(s1);
%     s1 = detrend( s1 ) ;
%     s1 = s1./std( s1 ) ;
%     m  = sqrt( mean( s1.^2, 2 ) )  ;      %%------------ spatial RMS
%     v  = var( s1, [], 2) ;              %%------------ spatial variance
%     mu = mean( s1, 2 ) ;                %%-----------spatia mean or ensemble mean
%     mu_t = mean( s1, 1 ) ;                %%---------temporal mean
%     vt  = var( s1, [], 1 ) ;             %%---------temporal variance
%     mean_temp_var = mean( vt ) ;
%     sync_mu(i) = var( mu )/mean_temp_var ;        %%---------spatial synchrony
%     rho = cov( s1 ) ; sync_rho(i) = mean( mean( rho ) ) ; %%---------spatial synchrony
%     fc  = corr(s1); fc=fc.*~eye(size(fc));sync_fc(i) = mean(mean((fc)));
%     cv1 = cov( log( m ), log( v ) ) ;  cv1=cv1(2, 1);
%     b_actual(i) = cv1/var( log( m ) ) ;
%     lm  = fitlm( log10( m ), log10( v ) ) ;  r2(i) = lm.Rsquared.Ordinary ;
%     ab2 = table2array( lm.Coefficients ) ; p_val(i) = ab2( 2, 4 ) ;
%     b_fit(i) = ab2( 2, 1 ) ;
%     cv = cov(m, v) ;  cv = cv(2, 1) ;
%     n  = size(s1,2);
%     b_analytical(i) = ( ((n-1)/n) * mean(m) * cv )/( ( mean(m)^2 + var(m) -var(mu) ) * var(m) );
% end


function [b_fit,p_val,r2,sync,b_analytical]=TL_nki(bold)
% % % bold : input data of size: time x location
sub=size(bold,1);  % total number of subjects
p_val=zeros(sub,1); r2=zeros(sub,1); b_fit=zeros(sub,1); b_analytical=zeros(sub,1); sync=zeros(sub,1); 

for i=1:sub
    s1 = [] ;
    s1 = bold{i} ;
%     s1( :, 109 ) = [] ;
    s1 = detrend( s1 ) ;
    [b_fit(i),p_val(i),r2(i),sync(i),b_analytical(i)]=slop_sync(s1);
end