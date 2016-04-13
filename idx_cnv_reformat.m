function CNV=idx_snp_cnv_reformat( fname )
    
idxCnv=[];
CNV.chr = [];
CNV.start = [];
CNV.end = [];
fid=fopen(fname);
A=textscan(fid,'%f64%f64%f64');
cnvChr = A{1};
cnvStart = A{2};
cnvEnd = A{3};

fprintf( 1, 'Number of known CNV %.d\n', length(cnvChr) );
%A=dlmread(fname,'\t',1,0);
for currChr=1:23
    idxA = find(cnvChr==currChr);
    currStart = cnvStart(idxA);
    currEnd = cnvEnd(idxA);

    % Filter redundant CNV
    [uStart,idxU,J] = unique(currStart);
    currStart = uStart;
    currEnd = currEnd(idxU);

    [uEnd,idxU,J] = unique(currEnd);
    currEnd = uEnd;
    currStart = currStart(idxU);

    CNV.chr = [ CNV.chr; repmat( currChr, length(currStart), 1 ) ];
    CNV.start = [ CNV.start; sort(currStart) ];
    CNV.end = [ CNV.end; sort(currEnd) ];
end
