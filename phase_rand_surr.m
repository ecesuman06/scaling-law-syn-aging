function surr_bold = phase_rand_surr(bold_se,area_randomized)
% Input: bold signal, bold_se(t,n) --> rows are t different time-points;  columns are n different regions; Vector size : bold_se-->time x regions, 
%         noise_intensity-->default 0.5
%         area_randomized--->1,  univariate time series randomization
% Output: surrogated bold signal. Vector size : time x area
number_of_surr_output = 1 ;
[nfrms,nts]     = size(bold_se);
% [nfrms,nts]     = size(s1);
% area_randomized=1;
% noise_intensity=0.5;
if rem(nfrms,2) == 0;    nfrms       = nfrms - 1 ;    bold_se     = bold_se(1:nfrms, :) ;  end
len_ser         = (nfrms-1)/2 ;
interv1         = 2 : len_ser+1 ;
interv2         = len_ser + 2 : nfrms ;

fft_s           = fft( bold_se ) ;
% fft_s           = fft(s1 ) ;
% for number_of_surr_output = 1:1
ph_rnd          = 2*pi*rand( len_ser, nts/area_randomized ) ;

ph_interv1      = repmat(exp( 0.5*1i*ph_rnd ), 1, area_randomized ) ;
ph_interv2      = conj( flipud( ph_interv1)) ;

fft_s_surr      = fft_s ;

fft_s_surr(interv1, :) = fft_s(interv1, :).*ph_interv1 ;
fft_s_surr(interv2, :) = fft_s(interv2, :).*ph_interv2 ;

surr_bold(:, :, number_of_surr_output)   = real( ifft( fft_s_surr ));  % output
% end
% surr_bold_m=mean(surr_bold,3);
% surr_bold=squeeze(surr_bold);
