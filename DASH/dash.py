import astrodash

from os import listdir
from os.path import isfile, join

filenames = list(map(lambda x: './data/'+x, [f for f in listdir('./data/') if isfile(join('./data/', f))]))

classification = astrodash.Classify(filenames, classifyHost=False, smooth=6, rlapScores=False)

bestFits, redshifts, bestTypes, rejectionLabels, reliableFlags, redshiftErrs = classification.list_best_matches(n=5, saveFilename='DASH_matches.txt')
classification.plot_with_gui(indexToPlot=0)