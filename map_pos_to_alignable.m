function [READN,READT,WINDOWN,WINDOWT] = map_pos_to_alignable( READN, READT, WINDOWN, WINDOWT, alignableDir, chrs )

try

% Convert to ALIGNABLE position
for ci=1:length(chrs)
    chr=chrs(ci);
    alignableFile = [ alignableDir '/ALIGNABLE_' alignableDir '_chr' num2str(chr) '.mat' ];
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

catch
    fprintf(1,'ERROR: Improper format for alignable genome file\n');
end
