function [SEG,S2,BKP,POS,R,pval,CHR]=segment_solexa_logratios( READN, READT, chrLengthFile, chrs, W, aN, aT, p_bkp, p_merge, flagSpeedup )

%chrLengthFile = '/xchip/cancergenome/cancergenome04/Derek/solexa/code/chromInfo_hg18.txt';   % ARGUMENT
chrLengthFile

%% Load chromosome lengths
fid = fopen(chrLengthFile);
I = textscan(fid,'%u%f64%f64%s');
%fclose(fid);
chrLength = I{2};
alignableLength = I{3};
clear I;
alignableGenome = sum(alignableLength);
fclose(fid);

CHR = [];
R = [];
LWT = [];
RWT = [];
POS = [];
POSA = [];
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
%    alignablePos = sort(READT.alignable(in_c_T));   % SORT alignable positions

    [localRatio,idxN,cNT,cTT,lwT,rwT]=calc_log_ratio( normalPos, tumorPos, aN, aT, W );

    CHR = [ CHR; repmat(c,length(localRatio),1) ];
    R = [ R; localRatio ];
    LWT = [ LWT; lwT ];
    RWT = [ RWT; rwT ];
    IDXN = [ IDXN; idxN ];
    POS = [ POS; tumorPos ];
%    POSA = [ POSA; alignablePos ];

    CUMULN = [ CUMULN; cNT ];
    CUMULT = [ CUMULT; cTT ];
end
fprintf(1,'\n');
toc

tic
%%---  Calculate p-values  ---%
pval=calc_pval_lognormal_approx(W,LWT,RWT,R,flagSpeedup);
fprintf(1,'\n5) Calculate tumor read statistics:  ');
toc

tic
%BKP=filter_pval_v2( W, pval, R, CHR, POS, POSA, IDXN, p_bkp );
BKP=filter_pval_v2( W, pval, R, CHR, POS, IDXN, p_bkp );
fprintf(1,'\n6) Filter tumor initial bkp candidates:  ');
toc

tic
%SEG=initialize_seg_logratios( chrLength, alignableLength, CHR, BKP, aN, aT, POS, POSA, CUMULN, CUMULT );
SEG=initialize_seg_logratios( chrLength, CHR, BKP, aN, aT, POS, CUMULN, CUMULT );
fprintf(1,'\n7) Initialize TUMOR segment statistics:  ');
toc

tic
%S2 = merge_segments_logratio( BKP, SEG, alignableLength, aN, aT, p_merge );
S2 = merge_segments_logratio( BKP, SEG, aN, aT, p_merge );
fprintf(1,'\n8) Merge tumor segments:  ');
toc

