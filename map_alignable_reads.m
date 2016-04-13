function [READN,READT,WINDOWN,WINDOWT,RATIOS] = map_alignable_reads( sampleinfo, tumorName, alignableDir, windowFile, minNormalCount, medianCenter, germlineCNV, removeCNV, calcRatios )
%------------------------------------------------------------------------%
%  FILE: map_alignable_reads.m                                           %
%------------------------------------------------------------------------%
%                  SOFTWARE COPYRIGHT NOTICE AGREEMENT                   %
%      This software and its documentation are copyright (2008) by the   %
%  Broad Institute/Massachusetts Institute of Technology.  All rights    %
%  are reserved.  This software is supplied without any warranty or      %
%  guaranteed support whatsoever. Neither the Broad Institute nor MIT    %
%  can be responsible for its use, misuse, or functionality.             %
%------------------------------------------------------------------------%
%  PARAMETER LIST                                                        %
%    alignableDir Directory of Matlab files with alignable coordinates   %
%                   (e.g., ALIGNABLE_hg18_N36_D2_chr1.mat)               %
%                   Note that file names must include directory name     %
%                   in the format ALIGNABLE_[dirName]_chr1.mat           %
%    windowFile   Matlab structure with start coordinates of alignable   %
%                   windows, e.g. HG18_N36_D2_WINDOWS_100K.mat was       %
%                   calculated with ARACHNE parameters N=36, D=2         %
%    sampleinfo   Sample information file                                %
%    tumorName    Unique name for tumor sample                           %
%    minNormalCount  Minimum number of normal reads in above windows     %
%                    (default=1)                                         %
%    medianCenter Flag to center copy-number ratios to 1 (default=0)     %
%    germlineCNV  File of germline copy-number variants                  %
%                   (default='CNV.verified_080606.combined_reformat.txt')%
%    removeCNV    Flag to ignore copy-number variants  (default=1)       %
%    calcRatios   Flag to calculate copy-number ratios in 100kb alignable%
%                   windows  (default=0)                                 %
%------------------------------------------------------------------------%
%  OUTPUT VARIABLES                                                      %
%    matfile         Output Matlab file                                  %
%    ratiofile       Output tab-delimited text file with 100kb ratios    %
%                                                                        %
%    READN,READT     Coordinates of aligned sequence reads               %
%      .chr          Chromosome                                          %
%      .pos          Start position of sequence alignment                %
%      .alignable    Start position in alignable coordinates             %
%                                                                        %
%    WINDOWN,WINDOWT Counts of aligned reads in 100kb alignable windows  %
%      .chr          Chromosome                                          %
%      .counts       Number of reads in 100kb alignable windows          %
%      .breaks       Start positions of 100 kb alignable windows         %
%      .alignable    Start positions in alignable coordinates            %
%                                                                        %
%    RATIOS          Tumor vs. normal copy-number ratios                 %
%      .chr          Chromosome                                          %
%      .windows      Start positions of 100 kb alignable windows         %
%      .windowIndex  Unique index of 100 kb alignable windows            %
%      .ratios       Copy-number ratio between tumor vs. normal          %
%------------------------------------------------------------------------%

format compact;

% Set names of output files
matfile = [ tumorName '_aligned_reads.mat' ]
ratiofile = [ tumorName '_alignable_100kb_ratios.txt' ];


%%---  LOAD DATA  ---%%
% WINDOWS of equal size in alignable portion of hg18 genome build
ALIGN=load(windowFile);

%% LOAD positions of copy number variants
if removeCNV
    CNV = idx_cnv_reformat(germlineCNV);
end

% Load info file with .marks.visual filenames in first column
fid = fopen(sampleinfo);
if fid == -1
    error(['ERROR: Cannot open sample information file ' sampleinfo]);
end
I = textscan(fid,'%s%s%s%*[^\n]','headerLines',1);  % First row is header information
fclose(fid);
FILES.file = I{1};
FILES.sample = I{2};
FILES.type = I{3};
clear I;

sampleNames = sort(unique(FILES.sample));

tumorIndex = find(strcmp(FILES.type,'Tumor')==1);
normalIndex = find(strcmp(FILES.type,'Normal')==1);
if length(tumorIndex) < 1
    error('ERROR: Second column must contain at least one Tumor');
end
if length(normalIndex) < 1
    error('ERROR: Second column must contain at least one Tumor');
end

% Load coordinates of aligned tumor reads
file=FILES.file(tumorIndex(1))
[WINDOWT,READT] = read_align_pos_file( file, ALIGN.WINDOWS );
READT.lane = repmat(1,length(READT.chr),1);

for i=2:length(tumorIndex)
    file = FILES.file(tumorIndex(i))
    [W,R] = read_align_pos_file( file, ALIGN.WINDOWS );

    WINDOWT.counts = WINDOWT.counts + W.counts;
    READT.chr = [ READT.chr; R.chr ];
    READT.pos = [ READT.pos; R.pos ];
    READT.lane = [ READT.lane; repmat(i,length(R.chr),1) ];
end

% Load coordinates of aligned normal reads
file = FILES.file(normalIndex(1))
[WINDOWN,READN] = read_align_pos_file( file, ALIGN.WINDOWS );
READN.lane = repmat(1,length(READN.chr),1);

for i=2:length(normalIndex)
    file = FILES.file(normalIndex(i))
    [W,R] = read_align_pos_file( file, ALIGN.WINDOWS );

    WINDOWN.counts = WINDOWN.counts + W.counts;
    READN.chr = [ READN.chr; R.chr ];
    READN.pos = [ READN.pos; R.pos ];
    READN.lane = [ READN.lane; repmat(i,length(R.chr),1) ]; 
end

%%---  SORT READ POSITIONS  ---%%
compound = READN.chr * 1e9 + READN.pos;
[c,readOrder] = sort(compound);
READN.chr = READN.chr(readOrder);
READN.pos = READN.pos(readOrder);
READN.lane= READN.lane(readOrder);
idxN = find(READN.chr>0);
READN.chr = READN.chr(idxN);
READN.pos = READN.pos(idxN);
READN.lane= READN.lane(idxN);

compound = READT.chr * 1e9 + READT.pos;
[c,readOrder] = sort(compound);
READT.chr = READT.chr(readOrder);
READT.pos = READT.pos(readOrder);
READT.lane= READT.lane(readOrder);
idxT = find(READT.chr>0);
READT.chr = READT.chr(idxT);
READT.pos = READT.pos(idxT);
READT.lane= READT.lane(idxT);

chrs=unique(READN.chr);

%%%--  MASK COPY NUMBER VARIANTS  ---%%
if removeCNV
  idxKeepN = ones(length(READN.chr),1);
  idxKeepT = ones(length(READT.chr),1);

  for ci=1:length(chrs)
    chr=chrs(ci);
    idxChrCNV = find(CNV.chr==chr);
    idxChrN = find(READN.chr==chr);
    idxChrT = find(READT.chr==chr);
    
    posN = READN.pos(idxChrN);
    posT = READT.pos(idxChrT);

    for i = 1:length(idxChrCNV)
	cnvL = CNV.start(idxChrCNV(i));
        cnvR = CNV.end(idxChrCNV(i));
        idxCurrN = find(posN>=cnvL & posN<=cnvR);
        idxCurrT = find(posT>=cnvL & posT<=cnvR);

        idxKeepN(idxChrN(idxCurrN)) = 0;
        idxKeepT(idxChrT(idxCurrT)) = 0;
    end
  end
    fprintf(1,'Normal reads in CNV:  %.0f\n', length(idxKeepN) - sum(idxKeepN) );
    fprintf(1,'Tumor  reads in CNV:  %.0f\n', length(idxKeepT) - sum(idxKeepT) );
end

READN.chr = READN.chr(find(idxKeepN==1));
READN.pos = READN.pos(find(idxKeepN==1));
READN.lane = READN.lane(find(idxKeepN==1));

READT.chr = READT.chr(find(idxKeepT==1));
READT.pos = READT.pos(find(idxKeepT==1));
READT.lane = READT.lane(find(idxKeepT==1));

%%---  MAP READ POSITIONS TO ALIGNABLE GENOME  ---%%
READN.alignable = [];
READT.alignable = [];
WINDOWN.alignable = [];
WINDOWT.alignable = [];

if ~calcRatios

    % Check availability of alignable genome information
    alignableFile = [ alignableDir '/ALIGNABLE_' alignableDir '_chr' num2str(chr) '.mat' ];
    try
        H = load(alignableFile);
        clear H;

        [READN,READT,WINDOWN,WINDOWT] = map_pos_to_alignable( READN, READT, WINDOWN, WINDOWT, alignableDir, chrs );

        RATIOS = make_ratios_struct( WINDOWN, WINDOWT, minNormalCount, medianCenter );

        % Write copy-number ratios file for 100 kb alignable windows
        fid=fopen( ratiofile, 'w' );
        fmt=['Chromosome\tStart\tCopy ratio\tTumor counts\tNormal counts\n'];
        fprintf(fid,fmt);
        fmt=['%d\t%d\t%.4f\t%d\t%d\n'];
        for i=1:length(RATIOS.ratios)
            fprintf(fid,fmt,RATIOS.chr(i),RATIOS.windows(i),RATIOS.ratios(i),WINDOWT.counts(i),WINDOWN.counts(i));
        end
        fclose(fid);

    % No alignable genome coordinates
    catch
        RATIOS = [];
        fprintf(1,['ERROR: No alignable genome coordinates found in directory ' alignableDir '\n']);
    end
else
    RATIOS = [];
end


%%---  SAVE OUTPUT  ---%
save(matfile,'READN','READT','RATIOS','WINDOWN','WINDOWT','sampleNames');

