function [ratios,diffs,pval]=update_seg_logratios( countN, countT, aN, aT, oldP, rmidx )
%------------------------------------------------------------------------%
%  FILE: estimate_seg_ratios.m                                           %
%------------------------------------------------------------------------%
%                  SOFTWARE COPYRIGHT NOTICE AGREEMENT                   %
%      This software and its documentation are copyright (2008) by the   %
%  Broad Institute/Massachusetts Institute of Technology.  All rights    %
%  are reserved.  This software is supplied without any warranty or      %
%  guaranteed support whatsoever. Neither the Broad Institute nor MIT    %
%  can be responsible for its use, misuse, or functionality.             %
%------------------------------------------------------------------------%
%  INPUT:  Current list of chromosomal breakpoints                       %
%  OUTPUT: Summary statistics for segments between these breakpoints     %
%------------------------------------------------------------------------%
%  PARAMETER LIST                                                        %
%    SEG_PARAM       Number of aligned reads in entire genome            %
%    chr             Chromosome number                                   %
%    bkpA            B+2 vector of breakpoints, plus chr start and end   %
%                      NOTE: Bkp are in ALIGNABLE genome coordinates     %
%    countN, countT  B+1 vector with counts of normal or tumor reads     %
%                      in segments flanking breakpoints                  %
%------------------------------------------------------------------------%
%  OUTPUT VARIABLES                                                      %
%    ratios          B+1 vector, copy-number ratios for the segments     %
%                                between adjacent breakpoints            %
%    diffs           B vector, differences of copy-number ratios         %
%                              between adjacent segments                 %
%    z               B vector, Z-scores between adjacent segments        %
%    lrt             B vector: sigma(left)*sigma(right)/sigma(both)      %
%    sd              B+1 vector for segments flanking breakpoints        %
%    spanSD          B vector, each SD assigned to a putative bkp        %
%    spanSizes       B vector, length of segment spanning bkp in         %
%                                    ALIGNABLE genome coordinates        %
%------------------------------------------------------------------------%

ratios = ( countT ./ aT ) ./ ( countN ./ aN );
diffs = log( ratios(2:end) ) - log( ratios(1:end-1) );

% NULL hypothesis: span corresponds to "spanning segments" without breakpoint
if length(countN) == 1                    % Only ONE spanning segment
    spanN = countN(1);
    spanT = countT(1);
else                                    
    spanN = countN(1:end-1) + countN(2:end);
    spanT = countT(1:end-1) + countT(2:end);
end

fracL = countN(1:end-1) ./ spanN;
fracR = countN(2:end) ./ spanN;
pval = ones(length(diffs),1);

if rmidx < 3 | rmidx > size(oldP) - 2 
    for b=1:length(diffs)
    [app_mu_left,app_sig_left(b)]=lognormal_approx_for_normratio( ...
    spanT(b)*fracL(b),sqrt(spanT(b)*fracL(b)),countN(b),sqrt(countN(b)));

    [app_mu_right,app_sig_right(b)]=lognormal_approx_for_normratio( ...
    spanT(b)*fracR(b),sqrt(spanT(b)*fracR(b)),countN(b+1),sqrt(countN(b+1)));
       tmp = normcdf(diffs(b),0,sqrt(app_sig_left(b)^2+app_sig_right(b)^2));
       pval(b) = 2*min(tmp,1-tmp);
   end
   pval = [ 1; pval ];
else
    b=rmidx-2;
    [app_mu_left,app_sig_left(b)]=lognormal_approx_for_normratio( ...
    spanT(b)*fracL(b),sqrt(spanT(b)*fracL(b)),countN(b),sqrt(countN(b)));

    [app_mu_right,app_sig_right(b)]=lognormal_approx_for_normratio( ...
    spanT(b)*fracR(b),sqrt(spanT(b)*fracR(b)),countN(b+1),sqrt(countN(b+1)));
    tmp = normcdf(diffs(b),0,sqrt(app_sig_left(b)^2+app_sig_right(b)^2));
    updateP1=2*min(tmp,1-tmp);

    b=rmidx-1;
    [app_mu_left,app_sig_left(b)]=lognormal_approx_for_normratio( ...
    spanT(b)*fracL(b),sqrt(spanT(b)*fracL(b)),countN(b),sqrt(countN(b)));

    [app_mu_right,app_sig_right(b)]=lognormal_approx_for_normratio( ...
    spanT(b)*fracR(b),sqrt(spanT(b)*fracR(b)),countN(b+1),sqrt(countN(b+1)));
   tmp = normcdf(diffs(b),0,sqrt(app_sig_left(b)^2+app_sig_right(b)^2));
updateP2=2*min(tmp,1-tmp);


    pval = [ oldP(1:rmidx-2); updateP1; updateP2; oldP(rmidx+2:end) ];

end

diffs = [ 0 ; diffs ];
