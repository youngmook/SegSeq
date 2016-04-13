function [BKP,cutoff_p]=filter_pval_numFP(W,pval,chrs,physicalPos,idxN,numFP)


BKP.idx = [];
BKP.chr = [];
BKP.delta = [];
BKP.pos = [];
up=[];
dn=[];
i=1;

%%---  Look for GAINS  ---%
if(0)
[mx,mi]=min(pval);
startIdx = max([ mi-W 1]);
endIdx = min([mi+W length(pval)]);
localIdx = startIdx:endIdx;
idxTies = find(pval(localIdx)==mx);
mi=floor(median(localIdx(idxTies)));
end

[orderP,orderIdx] = sort(pval);
% Only look through the top 10% of p-values
maxIdx = floor(length(orderIdx)/10);
idxTestPos = orderIdx(1:maxIdx);
maskPos = zeros(length(idxTestPos),1);

i=1;  % Reset candidate breakpoint number
b=0;  % Initialize scanning index

% Repeat until the largest local difference statistic falls below cutoff
while size(up,1) < numFP
  b = b + 1;

  if maskPos(b)==0

    mi=orderIdx(b);
    cutoff_p = pval(mi);

    if size(up,1) > 0
      minDistN = min(abs( repmat(idxN(mi),size(up,1),1) - up(:,5) ) );
    else
      minDistN = idxN(mi);
    end

    % Keep breakpoint if at least W normal reads away from the previous gain
    if minDistN > W
      up(i,1)=mi;
      up(i,2)=physicalPos(mi);
      up(i,3)=pval(mi);
%      up(i,4)=alignablePos(mi);
      up(i,5)=idxN(mi);
      up(i,6)=chrs(mi);
    end   

    i = i + 1;     % Increment breakpoint number

    startIdx = max([ mi-W 1]);
    endIdx = min([mi+W length(pval)]);

    % Mask local region 
    maskPos(find(idxTestPos>startIdx & idxTestPos<endIdx)) = 1;
  end
end

%%---  Gather all candidate breakpoints  ---%
if ~isempty(up)
    p=up;
    [tmp,si]=sort(p(:,1));
    p=p(si,:);
    BKP.idx = p(:,1);
    BKP.pos = p(:,2);
    BKP.delta = p(:,3);
%    BKP.alignablepos = p(:,4);
    BKP.chr = p(:,6);
    
%    [uVal,uIdxA] = unique(BKP.alignablepos,'first');
%    uIdx=uIdxA(find(BKP.delta(uIdxA)~=0));
%    BKP.chr=BKP.chr(uIdx);
%    BKP.idx=BKP.idx(uIdx);
%    BKP.pos=BKP.pos(uIdx);
%    BKP.delta=BKP.delta(uIdx);
%    BKP.alignablepos=BKP.alignablepos(uIdx);
end

