#!/opt/anaconda2/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 10:32:38 2016

@author: awm018
"""

from __future__ import print_function

def linegen(fobj):
    for line in fobj:
        yield line
    pass


def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    from itertools import izip_longest
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)


def main():

    from collections import Counter
    hist = Counter()
    hist.update(range(150, 1400))
    hist.subtract(range(150, 1400))

    for fname in ['small5.minisam',
                  'small6.minisam',
                  'small7.minisam',
                  'small8.minisam',
                  'small9.minisam',
                  'small10.minisam',
                  'medium5.minisam',
                  'medium6.minisam',
                  'medium7.minisam',
                  'medium8.minisam',
                  'medium9.minisam',
                  'medium10.minisam',
                  'large5.minisam',
                  ]:
        with open(fname, 'r') as f:
            pgen = grouper(f, 2, "")
            # distance between start positions
            #distgen = (abs(int(tail.split(',')[2]) - int(head.split(',')[2])) for head, tail in pgen) # mu 450.111745779 sigma 34.1305247875

            # gap width
#            distgen = (
#                min(
#                    abs(int(tail.split(',')[3]) - int(head.split(',')[2])),
#                    abs(int(tail.split(',')[2]) - int(head.split(',')[3]))
#                )
#                for head, tail in pgen) # mu 301.189271986 sigma 33.7575634303

            # center distance (doubled) # even are more common than odd, mu 900.22108274 sigma 67.8883766252
            # center distance # mu 450.017915846 sigma 33.9447897372
            distgen = (
                abs(
                    int(tail.split(',')[3]) + int(tail.split(',')[2])
                    -
                    (int(head.split(',')[3]) + int(head.split(',')[2]))
                ) / 2
                for head, tail in pgen)
            hist.update(distgen)
            pass
        pass
    print(hist)
    print(sorted(hist.keys()))

    from scipy.stats import norm
    nmu, nsigma = norm.fit(list(hist.elements()))
    print("mu", nmu, "sigma", nsigma)
#    from scipy.stats import  cauchy
#    cmu, csigma = cauchy.fit(list(hist.elements()))

    import matplotlib.pyplot as plt
    plt.figure()
    from numpy import array
    values = array(hist.values(), dtype=float) / sum(hist.values())
    plt.bar(hist.keys(), values, alpha=0.4, color='r', edgecolor='r', linewidth=0, width=1.0)

    pdf = norm.pdf(hist.keys(), loc=nmu, scale=nsigma)
    plt.plot(hist.keys(), pdf, 'b-')
#    pdf = cauchy.pdf(hist.keys(), loc=cmu, scale=csigma)
#    plt.plot(hist.keys(), pdf, 'g-')

    plt.show()

    pass

if __name__ == "__main__":
    main()
    pass
