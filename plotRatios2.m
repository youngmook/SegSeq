function xsize=plotRatios2( RATIOS, SEGS, SEG2 )
%  plotRatios.m
%  INPUT:  RATIOS structures (WC.chr, WC.windows, WC.counts)
%  OUTPUT:  Figure with ratio of binomials along genome coordinates
%
%  Derek Chiang
%  dchiang@broad.mit.edu
%

%selectWindows = find(isnan(RATIOS.ratios)==0);
somatic=find(RATIOS.chr<23);
selectWindows=intersect(somatic,find(RATIOS.ratios>0.01));
ratios = RATIOS.ratios(selectWindows);
chr = RATIOS.chr(selectWindows);
pos = RATIOS.windows(selectWindows);
%somatic=find(RATIOS.chr<23);
%medianCorrect = 1 - median(RATIOS.ratios(somatic))
%ratios = RATIOS.ratios + medianCorrect;

chrs=unique(RATIOS.chr);
windowPos = 1:length(selectWindows);

plot(windowPos(ratios>0),ratios(ratios>0),'ko','MarkerSize',1,'MarkerFaceColor','k');
axis([0 length(selectWindows) 0 max(ratios) + 0.2]);
%axis([0 length(selectWindows) 0 9]);
hold on;

plot( [1 length(ratios)], [1 1], 'k' );
currChrStart = 1;
chrStarts = [ 1 ];
for c=2:22;
    currChrStart = min(find(chr==c));  % This line changed
    chrStarts = [ chrStarts currChrStart ];
    plot( [currChrStart currChrStart], [0 max(ratios)+8], ':k','LineWidth',1 );
end
%chrStarts(23) = min(find(chr==23));

if exist( 'SEGS', 'var' )
    for c=1:22
%    for currChr=1:length(chrs)
%        c=chrs(currChr);
        segMean = SEGS.ratios(SEGS.chr==c);
        numWindows = histc( RATIOS.windows(RATIOS.chr==c), SEGS.left(SEGS.chr==c) );
        cumWindows = [ 0; cumsum(numWindows) ];
        cumWindows = cumWindows(1:(end-1));
        for i=1:(length(cumWindows)-1)
%            [ RATIOS.windows(chrStarts(c)+cumWindows(i))/1e6 RATIOS.windows(chrStarts(c)+cumWindows(i+1))/1e6 segMean(i) ]
           plot( [chrStarts(c)+cumWindows(i) chrStarts(c)+cumWindows(i+1) ], ...
                  repmat(segMean(i), 1, 2), 'g-', 'LineWidth', 2.5);

        end
        if c==22
            chrEnd = length(RATIOS.windows);
        else
            chrEnd = chrStarts(c+1);
        end
        plot( [chrStarts(c)+cumWindows(end) chrEnd], ...
              repmat(segMean(length(cumWindows)), 1, 2), 'g-', 'LineWidth', 2.5);
    end

end

if exist( 'SEG2', 'var' )
    for c=1:22
%    for currChr=1:length(chrs)
%        c=chrs(currChr);
        segMean = SEG2.ratios(SEG2.chr==c);
        numWindows = histc( RATIOS.windows(RATIOS.chr==c), SEG2.left(SEG2.chr==c) );
        cumWindows = [ 0; cumsum(numWindows) ];
        cumWindows = cumWindows(1:(end-1));
        for i=1:(length(cumWindows)-1)
%            [ RATIOS.windows(chrStarts(c)+cumWindows(i))/1e6 RATIOS.windows(chrStarts(c)+cumWindows(i+1))/1e6 segMean(i) ]
            plot( [chrStarts(c)+cumWindows(i) chrStarts(c)+cumWindows(i+1) ], ...
                  repmat(segMean(i), 1, 2), 'b-', 'LineWidth', 2.5);
        end
        if c==22
            chrEnd = length(RATIOS.windows);
        else
            chrEnd = chrStarts(c+1);
        end
        plot( [chrStarts(c)+cumWindows(end) chrEnd], ...
              repmat(segMean(length(cumWindows)), 1, 2), 'b-', 'LineWidth', 2.5);
    end

end

hold off;


set(gca, 'xtick', chrStarts);
set(gca,'xticklabel',1:22);
xlabel('Chromosome');
ylabel('Solexa copy ratios');

xsize=length(selectWindows);
