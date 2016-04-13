function [localRatio,idxN,cNT,cTT,lwT,rwT]=calc_log_ratio( posN, posT, aN, aT, W )

%------------------------------------------------------------------------%
%  FILE: calc_local_diff.m                                               %
%------------------------------------------------------------------------%
%                  SOFTWARE COPYRIGHT NOTICE AGREEMENT                   %
%      This software and its documentation are copyright (2008) by the   %
%  Broad Institute/Massachusetts Institute of Technology.  All rights    %
%  are reserved.  This software is supplied without any warranty or      %
%  guaranteed support whatsoever. Neither the Broad Institute nor MIT    %
%  can be responsible for its use, misuse, or functionality.             %
%------------------------------------------------------------------------%
%  INPUT:  List of aligned positions for normal reads and tumor reads    %
%  OUTPUT: Local difference statistic and associated summaries           %
%------------------------------------------------------------------------%
%  PARAMETER LIST                                                        %
%    posN      Positions of normal reads from current chromosome         %
%    posT      Positions of tumor reads from current chromosome          %
%    aN        Total number of aligned reads from the normal sample      %
%    aT        Total number of aligned reads from the tumor sample       %
%    W         Size of local windows, i.e. # of consecutive normal reads %
%------------------------------------------------------------------------%
%  OUTPUT VARIABLES                                                      %
%    localDiff Local difference statistic                                %
%    idxN      Index of the closest normal read to each tumor read       %
%    cNT(i)    Number of tumor reads from position 1 to normal read i    %
%    cTT(j)    Number of tumor reads from position 1 to tumor read j     %
%    lwT(j)    Number of tumor reads in the LEFT window of tumor read j  %
%    rwT(j)    Number of tumor reads in the RIGHT window of tumor read j %
%------------------------------------------------------------------------%

  % Ensure reads are sorted in ascending order
  posN = sort(posN);
  posT = sort(posT);

  %%---  Find the closest NORMAL reads  ---%
  n = length(posT);

  % lwT, rwT: Tumor counts in left, right windows of TUMOR reads
  % Initialize variables
  lwT=ones(n,1);
  rwT=ones(n,1);
  lwN=ones(n,1);
  rwN=ones(n,1);

  pos = [ 0; posT ];         % Add left boundary for histogram bins

  hNT = histc(posN,pos);     % Number of normal reads between tumor reads
  hTT = histc(posT,pos);     % Number of tumor reads between tumor reads

  % First index indicates reads between telomere and first read
  % Last index indicates exact matches to last read
  cNT= cumsum(hNT);   % Cumulative sum of tumor reads between normal reads
  cTT= cumsum(hTT);   % Increments by definition, unless multiple reads

  leftNidx = cNT(1:end-1);
  rightNidx = cNT(1:end-1) + 1;
  leftNidx(find(leftNidx<1)) = 1;
  rightNidx(find(rightNidx>length(posN))) = cNT(end-1);

  flankDist = [ abs( posT - posN(leftNidx)) abs( posT - posN(rightNidx) ) ]; 
  flankIdx = [ leftNidx rightNidx ];
  
  [closDist,choiceN] = min(flankDist,[],2);

  % Index of CLOSEST normal
  idxN = flankIdx(:,1);
  idxN(find(choiceN==2)) = flankIdx(find(choiceN==2),2);
  

  %%---  Count the number of tumor reads between normal reads  ---%
  pos = [0; posN ];   % Histogram bins now anchored by NORMAL reads

  hNN = histc(posN,pos);
  hTN = histc(posT,pos);

  cNN= cumsum(hNN);   % Increments by definition, unless multiple reads
  cTN= cumsum(hTN);


  % k is defined from W+1 to (n-W)
  % k is an index on the TUMOR READS
  startIdx = min(find(idxN>W));
  endIdx = max(find(idxN<length(cTN)-W));
  for k=startIdx:endIdx
    % Need to count the number of TUMOR READS
    lwT(k) = cTN(idxN(k)) - cTN(idxN(k)-W); % Tumor reads in LEFT window

    % Should be W by definition, except for multiple reads
    lwN(k) = cNN(idxN(k)) - cNN(idxN(k)-W); % Normal reads in LEFT window

    rwT(k) = cTN(idxN(k)+W) - cTN(idxN(k)); % Tumor reads in RIGHT window
    rwN(k) = cNN(idxN(k)+W) - cNN(idxN(k)); % Normal reads in RIGHT window
  end
  
  r=(lwT./lwN)*(aN/aT);                   % Copy-number ratios on LEFT side
  localDiff=(rwT./rwN-lwT./lwN)*(aN/aT);  % Local difference statistic
  localRatio = log( rwT ./ lwT );

  cNT=cNT(1:(end-1));
  cTT=cTT(1:(end-1));
  
