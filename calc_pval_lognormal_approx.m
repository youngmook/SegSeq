function pval=calc_pval_lognormal_approx( W, LWT, RWT, R, flagSpeedup )

numOutlier = 10;

LRsum = LWT+RWT;

[uLR,uLRi,uLRj]=unique(LRsum);

u_sig=zeros(length(uLR),1);
pval=ones(length(R),1);

for i=1:length(uLR)
    [app_mu,u_sig_approx(i)]=lognormal_approx_for_normratio(uLR(i)/2,sqrt(uLR(i)/2),W,sqrt(W));
end


[suLRj,suLRji]=sort(uLRj);
posdif=[0; find(diff(suLRj)); length(uLRj)];

for i=1:length(uLR)
  % Find entries with the same tumor read sum
  pos_with_sum=suLRji((posdif(i)+1):posdif(i+1));
  [uR,uRi,uRj]=unique(R(pos_with_sum));  
  
  cur_p=ones(length(uR),1); 

  if flagSpeedup
    currSD = u_sig_approx(i) * sqrt(2);
    for j=1:length(uR)
      if ( abs(uR(j)) > 2 * currSD )
        tmp = normcdf(uR(j),0,currSD);
        cur_p(j)= 2*min(tmp,1-tmp);
      end
    end
    pval(pos_with_sum)=cur_p(uRj);
%    disp([ i length(pos_with_sum) length(uR)]);

  else
    for j=1:length(uR)
      tmp = normcdf(uR(j),0,sqrt(2)*u_sig_approx(i));
      cur_p(j)= 2*min(tmp,1-tmp);
    end
    pval(pos_with_sum)=cur_p(uRj);
%    disp([ i length(pos_with_sum) length(uR)]);
   end
end

