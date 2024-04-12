function script_miditoolbox(filename)

m = readmidi(filename)

figure
plotdist(pcdist1(m))

m = setvalues(m,'chan',1) % channel 0 does not work
figure
plotdist(ivdist1(m))

figure
plotdist(pcdist2(m))

figure, hold on
plotmelcontour(m,0.25,'abs',':r.')
plotmelcontour(m,1,'abs','-bo')

keyname(kkkey(m))
figure
keystrengths = kkcc(m)
plotdist(keystrengths)

m = setvalues(m,'dur',1) % duration needs to be positive to work
figure
onsetdist(m,4,'fig')

bestmeter = meter(m)

figure
onsetacorr(m,4,'fig')

figure
plothierarchy(m,'sec')

figure
segmentgestalt(m,'fig');

figure
segmentprob(m,.6,'fig')

%%

mus.score(filename,'Group')
