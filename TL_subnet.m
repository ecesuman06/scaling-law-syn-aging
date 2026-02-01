function [ b_fit, p_val, r2, sync, b_analytical ] = TL_subnet(bold, subnet)
% age=cell2mat(tm(:,3)); 
% bold=tm(:,2);
sub=size(bold,1);  % total number of subjects
p_val=zeros(sub,1); r2=zeros(sub,1); b_fit=zeros(sub,1); b_analytical=zeros(sub,1); sync=zeros(sub,1); 
for i=1:sub
    s2 = []; 
    s2 = bold{i};
    s2(:,109)=[];
    s2 = detrend(s2);
    s1 = s2(:, subnet) ;
    [b_fit(i),p_val(i),r2(i),sync(i),b_analytical(i)]=slop_sync(s1);
end