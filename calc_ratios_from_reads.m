function RATIOS=calc_ratios_from_reads( READN, READT, WINDOWS )

aN=length(find(READN.chr>0 & READN.chr < 24 ));
aT=length(find(READT.chr>0 & READT.chr < 24 ));

RATIOS.chr = [];
RATIOS.windows = [];
RATIOS.ratios = [];
RATIOS.countN = [];
RATIOS.countT = [];

%%---  Count reads within each window  ---%
for c=1:23
   fprintf(1,[num2str(c) '..']);

   currWindows = WINDOWS.breaks(find(WINDOWS.chr==c));
   currCountN = histc(READN.pos(find(READN.chr==c)), currWindows );
   currCountT = histc(READT.pos(find(READT.chr==c)), currWindows );

   % Last bin of histc are exact matches to final edge
   currCountN = currCountN(1:(end-1));    
   currPos = currWindows(1:length(currCountN));

   normalObs = find(currCountN > 0);   % MINIMUM COUNTS in normals
   propN = currCountN(normalObs) ./ aN;
   propT = currCountT(normalObs) ./ aT;
   windows = currPos(normalObs);

   RATIOS.chr = [ RATIOS.chr; repmat(c,length(windows),1) ];
   RATIOS.windows = [ RATIOS.windows; windows ];
   RATIOS.ratios = [ RATIOS.ratios; propT ./ propN ];
   RATIOS.countN = [ RATIOS.countN; currCountN(normalObs) ];
   RATIOS.countT = [ RATIOS.countT; currCountT(normalObs) ];
end
fprintf(1,'\n');
