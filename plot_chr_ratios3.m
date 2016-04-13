function region = plot_chr_ratios3( RATIOS, chr, left, right, SEG, tumorName, color )
%  plotRatios.m
%  INPUT:  RATIOS structures (WC.chr, WC.windows, WC.counts)
%  OUTPUT:  Figure with ratio of binomials along genome coordinates
%
%  Derek Chiang
%  dchiang@broad.mit.edu
%

somatic=find(RATIOS.chr<23);
%medianCorrect = 1 - median(RATIOS.ratios(somatic))
%ratios = RATIOS.ratios + medianCorrect;
selectChr = find(RATIOS.chr==chr);


if ~exist('left')
   left = 0; 
end
if ~exist('right')
    right=max(RATIOS.windows(selectChr));
end
if ~exist('color')
    color='g';
end

if ~exist('tumorName')
    tumorName = '';
end

iRegion = intersect(find(RATIOS.windows(selectChr) >= left),...
                  find(RATIOS.windows(selectChr) <= right));
region = selectChr(iRegion);
              
%scatter(RATIOS.windows(region)/1e6,ratios(region),20,'bo','filled');
plot(RATIOS.windows(region)/1e6,RATIOS.ratios(region),'ko',...
    'MarkerSize',1.65,'MarkerFaceColor','k');
%axis([min(RATIOS.windows(region)/1e6)-0.2 max(RATIOS.windows(region)/1e6)+0.2 0 max(RATIOS.ratios(region))+0.2]);
%axis([min(RATIOS.windows(region)/1e6)-0.2 max(RATIOS.windows(region)/1e6)+0.2 0 max([2 ceil(max(RATIOS.ratios(region)))])]);
hold on;
%lot( [0 length(RATIOS.ratios)], [1 1], 'k' );
%for i=1:length(region)
%    plot( repmat(RATIOS.windows(region(i))/1e6,1,2), ...
%        [ ratios(region(i)) - 2*RATIOS.ci(region(i)) ...
%          ratios(region(i)) + 2*RATIOS.ci(region(i)) ],'k' );
%end

if exist( 'SEG', 'var' )
    idxSeg = find(SEG.chr==chr);
    chrLeft = SEG.left(idxSeg);
    chrRight = SEG.right(idxSeg);
    chrRatios = SEG.ratios(idxSeg);
    iRegion2 = intersect(find(SEG.left(idxSeg) >= left), ...
                find(SEG.right(idxSeg) <= right));
    for i=1:length(iRegion2)    
        plot( [chrLeft(iRegion2(i))/1e6 chrRight(iRegion2(i))/1e6], ...
           repmat(chrRatios(iRegion2(i)), 1, 2), [color '-'], 'LineWidth', 2);
        if i < length(iRegion2)
	    plot( repmat(chrRight(iRegion2(i)),1,2) / 1e6, ...
		  [ chrRatios(iRegion2(i)) chrRatios(iRegion2(i+1)) ], ...
		  [color '-'], 'LineWidth', 0.25 );
	end
    end
end

if exist( 'tumorName', 'var' )
    T=title([tumorName ' Chr ' num2str(chr)]);
    if chr==23
	T=title([tumorName ' Chr X']);
    end
    set(T,'FontName','Arial');
    set(T,'FontWeight','bold');
end

hold off;


X=xlabel(['Chr ' num2str(chr) ' position (Mb)']);
if chr==23
    X=xlabel(['Chr X position (Mb)']);
end
%xlabel(['Genomic position (Mb)']);
Y=ylabel('Copy-number ratio');
region = region;
%regionCI = RATIOS.ci(region);
set(X,'FontName','Arial');
set(Y,'FontName','Arial');
