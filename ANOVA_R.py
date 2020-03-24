import pandas as pd
import numpy as np
import scipy.stats as stats
from statsmodels.stats.libqsturng import psturng
import numpy as np

d = {}
d['a'] = [1, 1, 1,1,2,2,2,1,1,2,1,3,1,3,2,1,1,1,1,1]
d['aa'] = [1, 1, 1,1,2,2,2,1,1,2,1,3,1,3,2,1,1,1,1,1,1,1]
d['aaa'] = [1, 1, 1,1,2,2,2,1,1,2,1,3,1,3,2,1,1,1,1,1,2,2,2,1]
d['b'] = [3,3,3,2,2,2,3,3,3,3,2,2,2,3,4,8,8,8,8,8,8,12,12,12,12]
d['c'] = [2, 2, 2,2,2,2,2,2,2,2,1,1,1,2,2,2,2]
d['d'] = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
d['e'] = [2,2,2,2,2,2,2,2,2,2,2,2,5,6,4,6]
d['f'] = [1,2,3,4,5,6,7,8,9,1,2,3,4,5,6]
d['g'] = [1,1,1,1,1,1,1,1,1,1,10,1,1,1,1]
print(d)

subpopulations = ['a','b','c','d','e']
alphaSub = .005

def tukeyTtest(curDict, sub):
    td = {}
    if sub:
        subpop = [pop for pop in subpopulations]
        for pop1 in subpopulations:
            subpop.pop(0)
            for pop2 in subpop:
                T, P = stats.ttest_ind(np.asarray(curDict[pop1]), np.asarray(curDict[pop2]))
                if P < alphaSub:
                    td[P] = [pop1, pop2]
    return td


print(tukeyTtest(d, True))
