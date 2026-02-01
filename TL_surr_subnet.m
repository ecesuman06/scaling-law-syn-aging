function [b_fit,p_val,r2,sync,b_analytical]=TL_surr_subnet(bold,subnet)

sub=size(bold,1);
p_val=zeros(sub,1); r2=zeros(sub,1);b_fit=zeros(sub,1); b_analytical=zeros(sub,1); sync=zeros(sub,1); 

for i=1:sub
    s2 = []; 
    s2 = bold{i};
    s2(:,109)=[];
    s1 = s2(:, subnet) ;
    s1=phase_rand_surr(s1,1);
%     s1 = PhaseRand_surrogates(s1, 1);
    s1 = detrend(s1);
 [b_fit(i),p_val(i),r2(i),sync(i),b_analytical(i)]=slop_sync(s1);
end