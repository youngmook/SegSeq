function [BKP,POS,R,pval,p_bkp,p_merge]=segment_solexa_logratios_normals( READN, READT, chrLengthFile, chrs, W, aN, aT, numInitFP, numFinalFP, flagSpeedup )

%chrLengthFile = '/xchip/cancergenome/cancergenome04/Derek/solexa/code/chromInfo_hg18.txt';   % ARGUMENT

%% Load chromosome lengths
fid = fopen(chrLengthFile);
I = textscan(fid,'%u%f64%f64%s');
fclose(fid);
chrLength = I{2};
alignableLength = I{3};
clear I;


CHR = [];
R = [];
LWT = [];
RWT = [];
POS = [];
%POSA = [];
IDXN = [];
CUMULN = [];
CUMULT = [];

tic

for ci=1:length(chrs)
    c=chrs(ci);
    fprintf(1,[ num2str(c) '..' ]);
  
    in_c_N=find(READN.chr==c);
    in_c_T=find(READT.chr==c);
    normalPos = sort(READN.pos(in_c_N));            % SORT normal reads
    [tumorPos,tumorOrd] = sort(READT.pos(in_c_T));  % SORT tumor reads

    [localRatio,idxN,cNT,cTT,lwT,rwT]=calc_log_ratio( normalPos, tumorPos, aN, aT, W );

    CHR = [ CHR; repmat(c,length(localRatio),1) ];
    R = [ R; localRatio ];
    LWT = [ LWT; lwT ];
    RWT = [ RWT; rwT ];
    IDXN = [ IDXN; idxN ];
    POS = [ POS; tumorPos ];

    CUMULN = [ CUMULN; cNT ];
    CUMULT = [ CUMULT; cTT ];
end

%%---  Calculate p-values of local windows  ---%
LRsum=LWT+RWT;
pval=calc_pval_lognormal_approx(W,LWT,RWT,R,flagSpeedup);

if(0)
[uLR,uLRi,uLRj]=unique(LRsum);

u_sig=zeros(length(uLR),1);
pval=ones(length(R),1);
for i=1:length(uLR)
  [app_mu,u_sig(i)]=lognormal_approx_for_normratio(uLR(i)/2,sqrt(uLR(i)/2),W,sqrt(W));
  
  pos_with_sum=find(uLRj==i);
  [uR,uRi,uRj]=unique(R(pos_with_sum));  
  
  cur_p=ones(length(uR),1);

  for j=1:length(uR)
    currSD = u_sig(i) * sqrt(2);

    % Only calculate p-values when R is at least 2 std dev away from 0
      if ( abs(uR(j)) > 2 * currSD )
        tmp = normcdf(uR(j),0,currSD);
        cur_p(j)= 2*min(tmp,1-tmp);
      end

  end
  pval(pos_with_sum)=cur_p(uRj);
end


end
fprintf(1,'\n1) Normal read statistics:  ');
toc

tic
[BKP,p_bkp]=filter_pval_numFP( W, pval, CHR, POS, IDXN, numInitFP );
fprintf(1,'\n2) Filter 1000 FP initial normal breakpoints:  ');
toc

tic
  SEG=initialize_seg_logratios( chrLength, CHR, BKP, aN, aT, POS, CUMULN, CUMULT );
fprintf(1,'\n3) Initialize normal segments:  ');
toc

numBkp = 3e9;        % Dummy initialization variable
p_merge = p_bkp;

tic
while( numBkp > numFinalFP & p_merge > 1e-20 )
    S2 = merge_segments_logratio( BKP, SEG, aN, aT, p_merge );
    numBkp = length(S2.chr) - 23;

    p_merge = p_merge / 100
    SEG=S2;          % Start next iteration for p-value with segments so far
end

% Refine p-value
    p_merge = p_merge * 10;
    fprintf(1,[num2str(p_merge) '..'] );
    S2 = merge_segments_logratio( BKP, SEG, aN, aT, p_merge );
    numBkp = length(S2.chr) - 23;

    if ( numBkp > numFinalFP )
        p_merge = p_merge / 3;
    else
      p_merge = p_merge * 3;
    end

    fprintf(1,['\n4) Merge normal to ' num2str(numFinalFP) ' FP segments:  ']);
toc

