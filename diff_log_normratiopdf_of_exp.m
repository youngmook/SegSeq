function [d,d2]=diff_log_normratiopdf_of_exp(y,mu1,sig1,mu2,sig2)
p=pi;

expy=exp(y);
factor1=mu1.*expy./sig1.^2+mu2./sig2.^2;
factor2=expy.^2./sig1.^2+1./sig2.^2;
erf_factor=erf((factor1)./(factor2).^(1./2));
factor3=1./2.*mu1.^2./sig1.^2;
factor4=1./2.*mu2.^2./sig2.^2;
factor5=exp(1./2.*(factor1).^2./(factor2)-factor3-factor4);
factor6=exp(-factor3-factor4);
factor7=factor5./(factor2).^(3./2).*2.^(1./2)./p.^(1./2);

d=(expy.*(1./2.*(factor1).*factor7./sig1./sig2.*erf_factor+1./(factor2)./p./sig1./sig2.*factor6)+expy.*(1./2.*mu1.*expy./sig1.^3.*factor7./sig2.*erf_factor+1./2.*(factor1).*((factor1)./(factor2).*mu1.*expy./sig1.^2-(factor1).^2./(factor2).^2.*expy.^2./sig1.^2).*factor7./sig1./sig2.*erf_factor-3./2.*(factor1).*factor5./(factor2).^(5./2).*2.^(1./2)./p.^(1./2)./sig1.^3./sig2.*erf_factor.*expy.^2+(factor1).*factor7./sig1./sig2./pi.^(1./2).*exp(-(factor1).^2./(factor2)).*(mu1.*expy./sig1.^2./(factor2).^(1./2)-(factor1)./(factor2).^(3./2).*expy.^2./sig1.^2)-2./(factor2).^2./p./sig1.^3./sig2.*factor6.*expy.^2))./expy./(1./2.*(factor1).*factor7./sig1./sig2.*erf_factor+1./(factor2)./p./sig1./sig2.*factor6);

if nargout>1
d2=(expy.*(1./2.*(factor1).*exp(1./2.*(factor1).^2./ ...
                                                  (factor2)-factor3-factor4)./(factor2).^(3./2).*2.^(1./2)./p.^(1./2)./sig1./sig2.*erf_factor+1./(factor2)./p./sig1./sig2.*factor6)+2.*expy.*(1./2.*mu1.*expy./sig1.^3.*factor7./sig2.*erf_factor+1./2.*(factor1).*((factor1)./(factor2).*mu1.*expy./sig1.^2-(factor1).^2./(factor2).^2.*expy.^2./sig1.^2).*factor7./sig1./sig2.*erf_factor-3./2.*(factor1).*factor5./(factor2).^(5./2).*2.^(1./2)./p.^(1./2)./sig1.^3./sig2.*erf_factor.*expy.^2+(factor1).*factor7./sig1./sig2./pi.^(1./2).*exp(-(factor1).^2./(factor2)).*(mu1.*expy./sig1.^2./(factor2).^(1./2)-(factor1)./(factor2).^(3./2).*expy.^2./sig1.^2)-2./(factor2).^2./p./sig1.^3./sig2.*factor6.*expy.^2)+expy.*(1./2.*mu1.*expy./sig1.^3.*factor7./sig2.*erf_factor+mu1.*expy./sig1.^3.*((factor1)./(factor2).*mu1.*expy./sig1.^2-(factor1).^2./(factor2).^2.*expy.^2./sig1.^2).*factor7./sig2.*erf_factor-3.*mu1.*expy.^3./sig1.^5.*factor5./(factor2).^(5./2).*2.^(1./2)./p.^(1./2)./sig2.*erf_factor+2.*mu1.*expy./sig1.^3.*factor7./sig2./pi.^(1./2).*exp(-(factor1).^2./(factor2)).*(mu1.*expy./sig1.^2./(factor2).^(1./2)-(factor1)./(factor2).^(3./2).*expy.^2./sig1.^2)+1./2.*(factor1).*(mu1.^2.*expy.^2./sig1.^4./(factor2)-4.*(factor1)./(factor2).^2.*mu1.*expy.^3./sig1.^4+(factor1)./(factor2).*mu1.*expy./sig1.^2+4.*(factor1).^2./(factor2).^3.*expy.^4./sig1.^4-2.*(factor1).^2./(factor2).^2.*expy.^2./sig1.^2).*factor7./sig1./sig2.*erf((factor1)./(factor2).^(1./2))+1./2.*(factor1).*((factor1)./(factor2).*mu1.*expy./sig1.^2-(factor1).^2./(factor2).^2.*expy.^2./sig1.^2).^2.*factor7./sig1./sig2.*erf((factor1)./(factor2).^(1./2))-3.*(factor1).*((factor1)./(factor2).*mu1.*expy./sig1.^2-(factor1).^2./(factor2).^2.*expy.^2./sig1.^2).*factor5./(factor2).^(5./2).*2.^(1./2)./p.^(1./2)./sig1.^3./sig2.*erf((factor1)./(factor2).^(1./2)).*expy.^2+2.*(factor1).*((factor1)./(factor2).*mu1.*expy./sig1.^2-(factor1).^2./(factor2).^2.*expy.^2./sig1.^2).*factor7./sig1./sig2./pi.^(1./2).*exp(-(factor1).^2./(factor2)).*(mu1.*expy./sig1.^2./(factor2).^(1./2)-(factor1)./(factor2).^(3./2).*expy.^2./sig1.^2)+15./2.*(factor1).*factor5./(factor2).^(7./2).*2.^(1./2)./p.^(1./2)./sig1.^5./sig2.*erf((factor1)./(factor2).^(1./2)).*expy.^4-6.*(factor1).*factor5./(factor2).^(5./2).*2.^(1./2)./p.^(1./2)./sig1.^3./sig2./pi.^(1./2).*exp(-(factor1).^2./(factor2)).*(mu1.*expy./sig1.^2./(factor2).^(1./2)-(factor1)./(factor2).^(3./2).*expy.^2./sig1.^2).*expy.^2-3.*(factor1).*factor5./(factor2).^(5./2).*2.^(1./2)./p.^(1./2)./sig1.^3./sig2.*erf((factor1)./(factor2).^(1./2)).*expy.^2+(factor1).*factor7./sig1./sig2./pi.^(1./2).*(-2.*(factor1)./(factor2).*mu1.*expy./sig1.^2+2.*(factor1).^2./(factor2).^2.*expy.^2./sig1.^2).*exp(-(factor1).^2./(factor2)).*(mu1.*expy./sig1.^2./(factor2).^(1./2)-(factor1)./(factor2).^(3./2).*expy.^2./sig1.^2)+(factor1).*factor7./sig1./sig2./pi.^(1./2).*exp(-(factor1).^2./(factor2)).*(mu1.*expy./sig1.^2./(factor2).^(1./2)-2.*mu1.*expy.^3./sig1.^4./(factor2).^(3./2)+3.*(factor1)./(factor2).^(5./2).*expy.^4./sig1.^4-2.*(factor1)./(factor2).^(3./2).*expy.^2./sig1.^2)+8./(factor2).^3./p./sig1.^5./sig2.*factor6.*expy.^4-4./(factor2).^2./p./sig1.^3./sig2.*factor6.*expy.^2))./expy./(1./2.*(factor1).*factor7./sig1./sig2.*erf((factor1)./(factor2).^(1./2))+1./(factor2)./p./sig1./sig2.*factor6)-(expy.*(1./2.*(factor1).*factor7./sig1./sig2.*erf((factor1)./(factor2).^(1./2))+1./(factor2)./p./sig1./sig2.*factor6)+expy.*(1./2.*mu1.*expy./sig1.^3.*factor7./sig2.*erf((factor1)./(factor2).^(1./2))+1./2.*(factor1).*((factor1)./(factor2).*mu1.*expy./sig1.^2-(factor1).^2./(factor2).^2.*expy.^2./sig1.^2).*factor7./sig1./sig2.*erf((factor1)./(factor2).^(1./2))-3./2.*(factor1).*factor5./(factor2).^(5./2).*2.^(1./2)./p.^(1./2)./sig1.^3./sig2.*erf((factor1)./(factor2).^(1./2)).*expy.^2+(factor1).*factor7./sig1./sig2./pi.^(1./2).*exp(-(factor1).^2./(factor2)).*(mu1.*expy./sig1.^2./(factor2).^(1./2)-(factor1)./(factor2).^(3./2).*expy.^2./sig1.^2)-2./(factor2).^2./p./sig1.^3./sig2.*factor6.*expy.^2))./expy./(1./2.*(factor1).*factor7./sig1./sig2.*erf((factor1)./(factor2).^(1./2))+1./(factor2)./p./sig1./sig2.*factor6)-(expy.*(1./2.*(factor1).*factor7./sig1./sig2.*erf((factor1)./(factor2).^(1./2))+1./(factor2)./p./sig1./sig2.*factor6)+expy.*(1./2.*mu1.*expy./sig1.^3.*factor7./sig2.*erf((factor1)./(factor2).^(1./2))+1./2.*(factor1).*((factor1)./(factor2).*mu1.*expy./sig1.^2-(factor1).^2./(factor2).^2.*expy.^2./sig1.^2).*factor7./sig1./sig2.*erf((factor1)./(factor2).^(1./2))-3./2.*(factor1).*factor5./(factor2).^(5./2).*2.^(1./2)./p.^(1./2)./sig1.^3./sig2.*erf((factor1)./(factor2).^(1./2)).*expy.^2+(factor1).*factor7./sig1./sig2./pi.^(1./2).*exp(-(factor1).^2./(factor2)).*(mu1.*expy./sig1.^2./(factor2).^(1./2)-(factor1)./(factor2).^(3./2).*expy.^2./sig1.^2)-2./(factor2).^2./p./sig1.^3./sig2.*factor6.*expy.^2))./expy./(1./2.*(factor1).*factor7./sig1./sig2.*erf((factor1)./(factor2).^(1./2))+1./(factor2)./p./sig1./sig2.*factor6).^2.*(1./2.*mu1.*expy./sig1.^3.*factor7./sig2.*erf((factor1)./(factor2).^(1./2))+1./2.*(factor1).*((factor1)./(factor2).*mu1.*expy./sig1.^2-(factor1).^2./(factor2).^2.*expy.^2./sig1.^2).*factor7./sig1./sig2.*erf((factor1)./(factor2).^(1./2))-3./2.*(factor1).*factor5./(factor2).^(5./2).*2.^(1./2)./p.^(1./2)./sig1.^3./sig2.*erf((factor1)./(factor2).^(1./2)).*expy.^2+(factor1).*factor7./sig1./sig2./pi.^(1./2).*exp(-(factor1).^2./(factor2)).*(mu1.*expy./sig1.^2./(factor2).^(1./2)-(factor1)./(factor2).^(3./2).*expy.^2./sig1.^2)-2./(factor2).^2./p./sig1.^3./sig2.*factor6.*expy.^2);
 else
   d2=NaN;
end


