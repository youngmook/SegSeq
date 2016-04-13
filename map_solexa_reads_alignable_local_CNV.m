%function map_solexa_reads_alignable( alignableDir, outputDir, windowFile, sampleinfo, tumorName, minNormalCount, medianCenter )
%------------------------------------------------------------------------%
%  FILE: map_solexa_reads_alignable.m                                    %
%------------------------------------------------------------------------%
%                  SOFTWARE COPYRIGHT NOTICE AGREEMENT                   %
%      This software and its documentation are copyright (2008) by the   %
%  Broad Institute/Massachusetts Institute of Technology.  All rights    %
%  are reserved.  This software is supplied without any warranty or      %
%  guaranteed support whatsoever. Neither the Broad Institute nor MIT    %
%  can be responsible for its use, misuse, or functionality.             %
%------------------------------------------------------------------------%
%  PARAMETER LIST                                                        %
%    -a alignableDir Directory of Matlab files with alignable coordinates%
%                      (e.g., ALIGNABLE_hg18_N36_D2_chr1.mat)            %
%    -o outputDir    Directory for output Matlab & text files            %
%    -w windowFile   Matlab structure with start coordinates of alignable%
%                    windows, e.g. HG18_N36_D2_WINDOWS_100K.mat was      %
%                    calculated with ARACHNE parameters N=36, D=2        %
%    -i sampleinfo   Sample information file                             %
%    -n tumorName    Unique name for tumor sample                        %
%    -m minNormalCount  Minimum number of normal reads in above windows  %
%    -c medianCenter Flag to center copy-number ratios to 1 (default=0)  %
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


% PARAMETER DEFAULTS
useDefaults=1;
if useDefaults
    alignableDir = '/xchip/cancergenome04/Derek/solexa/hash_hg18/hg18_N36_D2';
    sampleinfodir = '/xchip/cancergenome04/Derek/solexa/sampleinfo/';
    outputDir = 'matfiles/'
    windowFile = 'HG18_N36_D2_WINDOWS_100K.mat';

%    sampleinfo = 'HCC1954_BL1954_solexa_sampleinfo.txt'
%    tumorName = 'HCC1954_noCNV'

%    sampleinfo = 'HCC1143_BL1143_5550_solexa_sampleinfo.txt'
%    tumorName = 'HCC1143_noCNV'

%    sampleinfo = 'H2347_BL2347_solexa_sampleinfo.txt';
%    tumorName = 'H2347_noCNV'

%      sampleinfo = 'BL1954_4343_2vs2_solexa_sampleinfo.txt'
%    tumorName = 'BL1954_2vs2_noCNV'

%    sampleinfo = 'BL1143_top2lanes_sampleinfo.txt'
%    tumorName = 'BL1143_top2lanes_noCNV'

      sampleinfo = 'BL2347_solexa_2vs2_sampleinfo.txt'
    tumorName = 'BL2347_2vs2_noCNV'


    minNormalCount = 1;
    medianCenter = 0;
end

matfile = [ outputDir tumorName '_aligned_reads.mat' ]
ratiofile = [ outputDir tumorName '_alignable_100kb_ratios.txt' ];

% WINDOWS of equal size in alignable portion of genome
ALIGN=load(windowFile);

%% LOAD positions of copy number variants
cnvHG18 = '/xchip/cancergenome04/Derek/solexa/power_calc/SNP6.0/CNV.verified_080606.combined_reformat.txt';
CNV = idx_cnv_reformat(cnvHG18);


% Load info file with .marks.visual filenames in first column
fid = fopen(sampleinfo);
I = textscan(fid,'%s%s%s%*[^\n]','headerLines',1);  % First row is header information
fclose(fid);
FILES.file = I{1};
FILES.sample = I{2};
FILES.type = I{3};
clear I;

sampleNames = sort(unique(FILES.sample));

% Load coordinates of aligned tumor reads
tumorIndex = find(strcmp(FILES.type,'Tumor')==1);
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
normalIndex = find(strcmp(FILES.type,'Normal')==1);
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

% SORT read positions
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


% MAP to alignable positions
READN.alignable = [];
READT.alignable = [];
WINDOWN.alignable = [];
WINDOWT.alignable = [];

numChr = length(unique(WINDOWN.chr));

%% MASK COPY NUMBER VARIANTS

idxKeepN = ones(length(READN.chr),1);
idxKeepT = ones(length(READT.chr),1);

for chr=1:numChr
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
%        keyboard

        idxKeepN(idxChrN(idxCurrN)) = 0;
        idxKeepT(idxChrT(idxCurrT)) = 0;
    end
end

    fprintf(1,'Normal reads in CNV:  %.0f\n', length(idxKeepN) - sum(idxKeepN) );
    fprintf(1,'Tumor  reads in CNV:  %.0f\n', length(idxKeepT) - sum(idxKeepT) );

READN.chr = READN.chr(find(idxKeepN==1));
READN.pos = READN.pos(find(idxKeepN==1));
READN.lane = READN.lane(find(idxKeepN==1));

READT.chr = READT.chr(find(idxKeepT==1));
READT.pos = READT.pos(find(idxKeepT==1));
READT.lane = READT.lane(find(idxKeepT==1));


% Convert to ALIGNABLE position
for chr=1:numChr

    alignableFile = [ alignableDir '/ALIGNABLE_chr' num2str(chr) '.mat' ];
    H = load(alignableFile);
    
    idxN = find(READN.chr==chr);
    phypos = READN.pos(idxN);
    phypos(find(phypos==0)) = 1;
    alignedN = H.ALIGNABLE.pos(phypos);
    
    idxT = find(READT.chr==chr);
    phypos = READT.pos(idxT);
    phypos(find(phypos==0)) = 1;
    alignedT = H.ALIGNABLE.pos(phypos);

    idxWN = find(WINDOWN.chr==chr);
    phypos = WINDOWN.breaks(idxWN);
    phypos(find(phypos==0)) = 1;
    alignedWN = H.ALIGNABLE.pos(phypos);

    idxWT = find(WINDOWT.chr==chr);
    phypos = WINDOWT.breaks(idxWT);
    phypos(find(phypos==0)) = 1;
    alignedWT = H.ALIGNABLE.pos(phypos);

    READN.alignable = [ READN.alignable; alignedN ];
    READT.alignable = [ READT.alignable; alignedT ];

    WINDOWN.alignable = [ WINDOWN.alignable; alignedWN ];
    WINDOWT.alignable = [ WINDOWT.alignable; alignedWT ];

    disp(['Chr ' num2str(chr) ': ' num2str(length(alignedN))]);
    clear H;
   
end

normalObs = find(WINDOWN.counts>=minNormalCount);   % MINIMUM COUNTS in normals
propT = WINDOWT.counts(normalObs) / sum(WINDOWT.counts(normalObs));
propN = WINDOWN.counts(normalObs) / sum(WINDOWN.counts(normalObs));
ratios = propT ./ propN;

somatic = find(WINDOWT.chr(normalObs)<23);

if medianCenter
    medianCorrect = 1 - median(ratios(somatic))
    ratios = ratios + medianCorrect;
end

RATIOS.propN = propN;
RATIOS.propT = propT;
RATIOS.ratios = ratios;
RATIOS.minNormalCount = minNormalCount;
RATIOS.windowIndex = normalObs;
RATIOS.windows = WINDOWT.breaks(normalObs);
RATIOS.chr = WINDOWT.chr(normalObs);
save(matfile,'READN','READT','RATIOS','WINDOWN','WINDOWT','sampleNames');

% Write 100kb alignable file
fid=fopen( ratiofile, 'w' );
fmt=['Chromosome\tStart\tCopy ratio\tTumor counts\tNormal counts\n'];
fprintf(fid,fmt);
fmt=['%d\t%d\t%.4f\t%d\t%d\n'];
for i=1:length(ratios)
     fprintf(fid,fmt,RATIOS.chr(i),RATIOS.windows(i),RATIOS.ratios(i),WINDOWT.counts(i),WINDOWN.counts(i));
end
fclose(fid);



