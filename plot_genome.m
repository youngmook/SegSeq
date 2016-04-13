function plot_genome(RATIOS,SEG,chrs,tumorName)
    for currChr=1:length(chrs)
        plotNum = mod(currChr,12);
        if plotNum == 1
            fh=figure(floor(currChr/4)+1);
            clf
	    set(fh,'Position',[0 0 800 1000]);
        end
        if plotNum == 0
            plotNum = 12;
        end
        subplot(6,2,plotNum);
plot_chr_ratios3(RATIOS,chrs(currChr),1,250e6,SEG,tumorName,'g');
    end
end
