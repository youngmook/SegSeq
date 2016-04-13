function write_segments_logratio( fname, SEG )

fid = fopen( fname, 'w' );

fmt=['Chromosome\tStart\tEnd\tCopy ratio\tp-value\tNormal count\tTumor count\n'];
fprintf(fid,fmt);
fmt=['%d\t%.0f\t%.0f\t%.4f\t%.2e\t%d\t%d\n'];

for i=1:length(SEG.chr)
    fprintf(fid, fmt, SEG.chr(i), SEG.left(i), SEG.right(i)-1, SEG.ratios(i), ...
            SEG.pval(i), SEG.countN(i), SEG.countT(i) );
end
fclose(fid);
