function WINDOWS = create_windows_from_reads( READN, W, chrBoundFile )
%
%  FILE: create_windows_from_uniform.m

chrBounds = dlmread(chrBoundFile);

WINDOWS.chr = [];
WINDOWS.breaks = [];

chrs=unique(READN.chr);
chrs=setdiff(chrs,0);
for ci=1:length(chrs)
    c=chrs(ci);

    in_chr=find(READN.chr==c);
    chrPos = READN.pos(in_chr);
    startPos = chrBounds(c,2);
    endPos = chrBounds(c,3);

    cenLpos = chrBounds(c,4);
    cenRpos = chrBounds(c,5);

    startIdx = min(find(chrPos>startPos));
    cenLidx = max(find(chrPos<cenLpos));
    cenRidx = min(find(chrPos>cenRpos));
    endIdx = max(find(chrPos<endPos));

    [ c startIdx cenLidx cenRidx endIdx ]

    currWindows = [ chrBounds(c,2); chrPos(startIdx:W:cenLidx); chrPos(cenRidx:W:endIdx); chrBounds(c,3) ];

    WINDOWS.breaks = [ WINDOWS.breaks; currWindows ];
    WINDOWS.chr = [ WINDOWS.chr; ones(length(currWindows),1,'single')*c ];
end
