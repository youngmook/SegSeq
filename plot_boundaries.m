function [normalReads,tumorReads]=plot_boundaries( READN, READT, currChr, searchStart, searchEnd, tumorName, color )


if ~exist( 'color', 'var' )
  color = 'c';
end

normalChrReads = READN.pos(find(READN.chr==currChr));
tumorChrReads = READT.pos(find(READT.chr==currChr));
    
normalReads = normalChrReads(intersect(find(normalChrReads>searchStart*1e6),...
              find(normalChrReads<searchEnd*1e6)));
tumorReads = tumorChrReads(intersect(find(tumorChrReads>searchStart*1e6),...
              find(tumorChrReads<searchEnd*1e6)));
normalReads = sort(normalReads);
tumorReads = sort(tumorReads);

%clf;
plot(normalReads/1e6, (1:length(normalReads))/length(normalReads), 'ko','MarkerFaceColor','k','MarkerSize',2);
hold on;
plot(tumorReads/1e6, (1:length(tumorReads))/length(tumorReads), [color 'o'],'MarkerFaceColor',color,'MarkerSize',2);
hold off;
axis([searchStart searchEnd 0 1]);
xlabel(['Chromosome ' num2str(currChr) ' position (Mb)']);
ylabel('Fraction of aligned reads in window');
title(tumorName);
legend('Normal','Tumor','Location','SouthEast');
