function RATIOS=make_ratios_struct( WINDOWN, WINDOWT, minNormalCount, medianCenter )

normalObs = find(WINDOWN.counts>=minNormalCount);   % MINIMUM COUNTS in normals
propT = WINDOWT.counts(normalObs) / sum(WINDOWT.counts(normalObs));
propN = WINDOWN.counts(normalObs) / sum(WINDOWN.counts(normalObs));
ratios = propT ./ propN;

somatic = find(WINDOWT.chr(normalObs)<23);

if medianCenter
    medianCorrect = 1 - median(ratios(somatic))
    ratios = ratios + medianCorrect;
end

RATIOS.propN = propN;
RATIOS.propT = propT;
RATIOS.ratios = ratios;
RATIOS.minNormalCount = minNormalCount;
RATIOS.windowIndex = normalObs;
RATIOS.windows = WINDOWT.breaks(normalObs);
RATIOS.chr = WINDOWT.chr(normalObs);
