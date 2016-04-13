function [S,removed_z]=merge_segments_logratio( BKP, SEG, aN, aT, p_merge )
%------------------------------------------------------------------------%
%  FILE: merge_segments.m                                                %
%------------------------------------------------------------------------%
%                  SOFTWARE COPYRIGHT NOTICE AGREEMENT                   %
%      This software and its documentation are copyright (2008) by the   %
%  Broad Institute/Massachusetts Institute of Technology.  All rights    %
%  are reserved.  This software is supplied without any warranty or      %
%  guaranteed support whatsoever. Neither the Broad Institute nor MIT    %
%  can be responsible for its use, misuse, or functionality.             %
%------------------------------------------------------------------------%
%  INPUT:  Initial list of copy-number segments and statistics           %
%          <- initialize_seg_ratios.m                                    %
%  OUTPUT: Merged list of copy-number segments with p-values < cutoff_p  %
%------------------------------------------------------------------------%
%  PARAMETER LIST                                                        %
%    BKP             Matlab structure with breakpoint information        %
%    SEG             Matlab structure with copy-number segment statistics%
%    SEG_PARAM       Matlab structure with cutoff values                 %
%    cutoff_p        Stopping point for merging adjacent segments        %
%    alignableLength Length of chromosomes, in alignable coordinates     %
%------------------------------------------------------------------------%
%  OUTPUT VARIABLES                                                      %
%    S               Matlab structure with copy-number segment statistics%
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


S.chr = [];
S.left = [];
S.right = [];
S.countN = [];
S.countT = [];
S.ratios = [];
S.diff = [];
S.pval = [];
%S.spanSizes = [];
removed_z = [];

chrs = unique(SEG.chr);
for ci=1:length(chrs)
    c = chrs(ci);
    
    idx = find(SEG.chr==c);
    idxChr = find(BKP.chr==c);

    currSeg = [ SEG.left(idx) SEG.right(idx) SEG.ratios(idx) SEG.diff(idx) SEG.pval(idx) SEG.countN(idx) SEG.countT(idx) ];

    countN = SEG.countN(idx);
    countT = SEG.countT(idx);
    
    if ( length(idx) > 1 )
      % Do not consider first segment in merging
      % Merge by p-values from least significant to most significant
      [maxP,mi] = max( currSeg(2:size(currSeg,1),5) );   % Z-score
      rmidx = mi+1;                      % Compensate for ignoring first segment

      % Continue merging segments until the lowest Z-score exceeds cutoff
      % or only one segment remains on the chromosome
      while maxP > p_merge & size(currSeg,1) > 1

                                       % Debug print statements
%        disp(['Remove bkp ' num2str(rmidx) ' chr ' num2str(c) ':' num2str(currSeg(rmidx,1)) '  p=' num2str(currSeg(rmidx,5)) '  N=' num2str(size(currSeg,1))]);


        removed_z = [ removed_z; currSeg(rmidx,5) ];
                                       % Update breakpoint coordinates
        if rmidx < size(currSeg,1)     % Merge an internal segment
            new_left = [ currSeg(1:(rmidx-1),1); currSeg(((rmidx+1):end),1)];
            new_right = [ currSeg(1:(rmidx-2),2); currSeg((rmidx:end),2)];

        else                           % Merge the last segment
            new_left = currSeg(1:(rmidx-1),1);
            new_right = [ currSeg(1:(rmidx-2),2); currSeg(rmidx,2) ];
        end

        % Remove LEFT breakpoint, so merge counts with previous segment
        countN = [ countN(1:(rmidx-2)); countN(rmidx-1) + countN(rmidx); countN((rmidx+1):end) ];
        countT = [ countT(1:(rmidx-2)); countT(rmidx-1) + countT(rmidx); countT((rmidx+1):end) ];
        
        emptyN = find(countN==0);
        for e=1:length(emptyN)
             disp(['Zero normal counts in chr ' num2str(c) ' at ' num2str(currSeg(e,1)) ' to ' num2str(currSeg(e,2))]);
        end

                                        % Calculate statistics for new segments
%        [ratios,diffs,pval] = estimate_seg_logratios( SEG_PARAM, countN, countT, sum(alignableLength) );
%        [ratios,diffs,pval] = update_seg_logratios( countN, countT, sum(alignableLength), aN, aT, currSeg(:,5), rmidx );
        [ratios,diffs,pval] = update_seg_logratios( countN, countT, aN, aT, currSeg(:,5), rmidx );


%	currSeg(rmidx-3:rmidx+3,5)'
%	pval(rmidx-3:rmidx+3)'

                                        % Update segment statistics
%        currSeg = [ new_left new_right ratios diffs pval spanSizes countN countT ];
        currSeg = [ new_left new_right ratios diffs pval countN countT ];

                                        % Repeat for next largest p-value
        if( size(currSeg,1) > 1 )
            [maxP,mi] = max( currSeg(2:size(currSeg,1),5) );  % p-value
            rmidx = mi+1;
        end
    end

    end

                                        % FINAL segments after merging
    S.chr = [S.chr; repmat(c,size(currSeg,1),1) ];
    S.left = [ S.left; currSeg(:,1) ];
    S.right = [ S.right; currSeg(:,2) ];
    S.ratios = [ S.ratios; currSeg(:,3) ];
    S.diff = [ S.diff; currSeg(:,4) ];
    S.pval = [ S.pval; currSeg(:,5) ];
    S.countN = [ S.countN; countN ];
    S.countT = [ S.countT; countT ];
end

