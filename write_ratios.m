function write_ratios( fname, RATIOS )

fid = fopen( fname, 'w' );

fmt=['Chromosome\tStart\tNormal counts\tTumor counts\tCopy ratio\n'];
fprintf(fid,fmt);
fmt=['%d\t%.0f\t%.0f\t%.0f\t%.3f\n'];

for i=1:length(RATIOS.chr)
    fprintf(fid, fmt, RATIOS.chr(i), RATIOS.windows(i), ...
	    RATIOS.countN(i), RATIOS.countT(i), RATIOS.ratios(i) );
end
fclose(fid);
