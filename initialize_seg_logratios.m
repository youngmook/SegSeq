function SEG=initialize_seg_logratios( chrLength, CHR, BKP, aN, aT, POS, CUMULN, CUMULT )
%------------------------------------------------------------------------%
%  FILE: initialize_seg_ratios.m                                         %
%------------------------------------------------------------------------%
%                  SOFTWARE COPYRIGHT NOTICE AGREEMENT                   %
%      This software and its documentation are copyright (2008) by the   %
%  Broad Institute/Massachusetts Institute of Technology.  All rights    %
%  are reserved.  This software is supplied without any warranty or      %
%  guaranteed support whatsoever. Neither the Broad Institute nor MIT    %
%  can be responsible for its use, misuse, or functionality.             %
%------------------------------------------------------------------------%
%  INPUT:  Initial list of candidate breakpoints for segmentation        %
%          <- filter_local_diff.m                                        %
%  OUTPUT: Initial list of segments and statistics                       %
%          -> merge_segments.m                                           %
%------------------------------------------------------------------------%
%  PARAMETER LIST                                                        %
%    chrLength       Length of chromosomes                               %
%    chrs            List of chromosomes for segmentation                %
%    BKP             Matlab structure with breakpoint information        %
%    SEG_PARAM       Matlab structure with cutoff values                 %
%    POS             Matlab structure with tumor read positions          %
%    CUMULN          Matlab structure with cumulative num of normal reads%
%    CUMULT          Matlab structure with cumulative num of tumor reads %
%------------------------------------------------------------------------%
%  OUTPUT VARIABLES                                                      %
%    SEG             Matlab structure with copy-number segment statistics%
%                      calculated from estimate_seg_ratios.m             %
%      .chr          Chromosome                                          %
%      .left         Segment start position (physical coordinates)       %
%      .right        Segment end position   (physical coordinates)       %
%      .countN       Number of aligned normal reads in segment           %
%      .countT       Number of aligned tumor reads in segment            %
%      .ratios       Copy-number ratios                                  %
%      .lrt          Likelihood ratio test statistic                     %
%      .sd           Gaussian approximation for standard deviation of    %
%                      copy-number ratio (current segment)               %
%      .diff         Local difference statistic between adjacent segments%
%      .z            Z-score approximation for copy-number differences   %
%      .spanSD       Gaussian approximation for standard deviation of    %
%                      copy-number ratio (two adjacent segments)         %
%      .spanSizes    Length of spanning segments (two adjacent segments) %
%------------------------------------------------------------------------%

SEG.chr = [];
SEG.left = [];
SEG.right = [];
SEG.countN = [];
SEG.countT = [];
SEG.ratios = [];
SEG.pval = [];
SEG.diff = [];
%SEG.spanSizes = [];

chrs=unique(CHR);
for ci=1:length(chrs)
    c=chrs(ci);

     idxChr = find(BKP.chr==c);
     maxIdx = max(find(CHR==c));

    % bkps vector includes the chromosome ends
    bkps = [ 1; BKP.pos(idxChr); chrLength(c) ];
    ends = [ POS(BKP.idx(idxChr) - 1); chrLength(c) ];
    
    countN = CUMULN( [ BKP.idx(idxChr); maxIdx ] ) - CUMULN( [ 1; BKP.idx(idxChr) ] );
    countT = CUMULT( [ BKP.idx(idxChr); maxIdx ] ) - CUMULT( [ 1; BKP.idx(idxChr) ] );
    
    % Segment statistics are evaluated in ALIGNABLE genome coordinates
    [ratios,diffs,pval] = estimate_seg_logratios( countN, countT, aN, aT );
    
    SEG.chr = [ SEG.chr; repmat(c, length(ratios), 1 ) ];
    SEG.left = [ SEG.left; bkps(1:(length(bkps)-1)) ];
    SEG.right = [ SEG.right; ends ];
    SEG.countN = [ SEG.countN; countN ];
    SEG.countT = [ SEG.countT; countT ];
    SEG.ratios = [ SEG.ratios; ratios ];
    SEG.pval = [ SEG.pval; pval ];
    SEG.diff = [ SEG.diff; diffs ];
%    SEG.spanSizes = [ SEG.spanSizes; spanSizes ];
end
