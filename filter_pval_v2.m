function BKP=filter_pval_v2(W,pval,diffR,chrs,physicalPos,idxN,cutoff_p)


BKP.idx = [];
BKP.chr = [];
BKP.delta = [];
BKP.pos = [];
up=[];
dn=[];
i=1;

pval_dn = pval;
pval_up = pval;

pval_dn(find(diffR>=0))=1;
pval_up(find(diffR<0))=1;

%profile on
%%---  Look for LOSSES  ---%
[orderP,orderIdx] = sort(pval_dn);
idxTestLoss = orderIdx(find(orderP <= cutoff_p));
maskLoss = zeros(length(idxTestLoss),1);


% Repeat until the smallest local difference statistic rises above cutoff
for b=1:length(maskLoss)
  if maskLoss(b) == 0
     mi = orderIdx(b);
     if size(dn,1) > 0 
         minDistN = min(abs( repmat(idxN(mi),size(dn,1),1) - dn(:,5) ) );
     else
         minDistN = idxN(mi);
     end

  % Keep breakpoint if at least W normal reads away from the previous loss
  if minDistN > W
    dn(i,1)=mi;
    dn(i,2)=physicalPos(mi);
    dn(i,3)=pval_dn(mi);
    dn(i,5)=idxN(mi); 
    dn(i,6)=chrs(mi);
  end

  i = i + 1;     % Increment breakpoint number
  % Remove local segment of width (2+W)
%  disp(['L Removing ' num2str(max([ mi-W 1])) ' to ' num2str(min([mi+W length(pval_dn)]))]);

   startIdx = max([ mi-W 1]);
   endIdx = min([mi+W length(pval_dn)]);

   % Mask local region 
   maskLoss(find(idxTestLoss>startIdx & idxTestLoss<endIdx)) = 1;
   end
end

%%---  Look for GAINS  ---%
[orderP,orderIdx] = sort(pval_up);
idxTestGain = orderIdx(find(orderP <= cutoff_p));
maskGain = zeros(length(idxTestGain),1);

i=1;  % Reset candidate breakpoint number

% Repeat until the largest local difference statistic falls below cutoff
for b=1:length(maskGain)
  if maskGain(b) == 0
     mi = orderIdx(b);
     if size(up,1) > 0 
         minDistN = min(abs( repmat(idxN(mi),size(up,1),1) - up(:,5) ) );
     else
         minDistN = idxN(mi);
     end

  % Keep breakpoint if at least W normal reads away from the previous gain
  if minDistN > W
    up(i,1)=mi;
    up(i,2)=physicalPos(mi);
    up(i,3)=pval_up(mi);
    up(i,5)=idxN(mi);
    up(i,6)=chrs(mi);
  end   

  i = i + 1;     % Increment breakpoint number
  % Remove local segment of width (2+W)
%  disp(['G Removing ' num2str(max([ mi-W 1])) ' to ' num2str(min([mi+W length(pval_dn)]))]);

   startIdx = max([ mi-W 1]);
   endIdx = min([mi+W length(pval_dn)]);

   % Mask local region 
   maskGain(find(idxTestGain>startIdx & idxTestGain<endIdx)) = 1;
   end
end

%profile report

%%---  Gather all candidate breakpoints  ---%
if ~isempty(up)
    p=[up;dn];
    [tmp,si]=sort(p(:,1));
    p=p(si,:);
    BKP.idx = p(:,1);
    BKP.pos = p(:,2);
    BKP.delta = p(:,3);
    BKP.idxN = p(:,5);
    BKP.chr = p(:,6);
end

