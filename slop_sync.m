    function [b_fit,p_val,R2,sync,b_analytical]=slop_sync(s1)

    s1 = s1./std( s1 ) ;

    m = sqrt(mean(s1.^2,2))  ;      %%------------ spatial RMS
    v  = var(s1,[],2);              %%------------ spatial variance

    mu = mean(s1,2);                %%-----------spatia mean or ensemble mean
%     mu_t=mean(s1,1);                %%---------temporal mean
    vt  = var(s1,[],1);             %%---------temporal variance
    mean_temp_var=mean(vt) ;        %%---------mean of temporal variance
    sync_mu=var(mu)/mean_temp_var;     sync=  sync_mu;            %%---------spatial synchrony
%     rho=cov(s1); sync_rho=mean(mean(rho)); sync=  sync_rho;     %%---------spatial synchrony
%     fc=corr(s1); sync_fc=mean(mean(abs(fc))); sync=  sync_fc;   %%---------spatial synchrony
    cv1 = cov(log(m), log(v)) ;  cv1=cv1(2, 1);
    b_actual = cv1/var(log(m));
    lm = fitlm(log10(m), log10(v)) ;  R2=lm.Rsquared.Ordinary ;
    ab2 = table2array(lm.Coefficients) ; p_val = ab2(2, 4) ;
    b_fit = ab2(2, 1);
    cv = cov(m, v) ;  cv = cv(2, 1) ;
    n=size(s1,2);
    b_analytical = ( ((n-1)/n) * mean(m) * cv )/( ( mean(m)^2 + var(m) -var(mu) ) * var(m) );