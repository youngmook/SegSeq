function xsize=plotRatios( RATIOS, gamma, mu1 )
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
%somatic=find(RATIOS.chr<23);
%medianCorrect = 1 - median(RATIOS.ratios(somatic))
%ratios = RATIOS.ratios + medianCorrect;


windowPos = 1:length(selectWindows);

plot(windowPos(ratios>0),ratios(ratios>0),'ko','MarkerSize',1,'MarkerFaceColor','k');
%axis([0 length(selectWindows) 0 max(ratios) + 0.2]);
axis([0 length(selectWindows) 0 2]);
hold on;

plot( [1 length(ratios)], [1 1], 'k' );
currChrStart = 1;
chrStarts = [ 1 ];
for c=2:22
    currChrStart = min(find(chr==c));  % This line changed
    chrStarts = [ chrStarts currChrStart ];
    plot( [currChrStart currChrStart], [0 max(ratios)+8], ':k','LineWidth',1 );
end
if exist( 'gamma', 'var' )
    [Y,iMax] = max(gamma,[],2);
    plot(1:length(selectWindows),mu1(iMax(selectWindows)),'b-','LineWidth',1.7);
end
hold off;


set(gca, 'xtick', chrStarts);
set(gca,'xticklabel',1:22);
xlabel('Chromosome');
ylabel('Solexa copy ratios');

xsize=length(selectWindows);