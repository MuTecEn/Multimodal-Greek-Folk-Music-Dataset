Examples:

Single file:
analyzeMIDI('Syrtos-Salamina.mid',[0,NaN,NaN,NaN,6],4)
(3rnd parameter is the number of beats per bar)

analyzeMIDI code detects when the file name is  Syrtos_Nysyros and performs a dedicated pitch spelling for that one.

Folder:
batchanalysis
You need to use my updated keysignature.csv file enclosed.

And then the tsne related to the pattern analysis:
patternembedding

Because the text display of the results is very long, exceeding the Command Windows limited buffer, you can print again the text display (and interrupt using ctrl+C) by using: textdisplay
