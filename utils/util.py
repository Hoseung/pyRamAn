# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 00:26:50 2015

@author: hoseung
"""
import os
import numpy as np
from glob import glob

def mkdir(dirpath):
    if not os.path.isdir(dirpath):
        os.mkdir(dirpath)

def point2range(cx,r):
    return [[cx[0]-r,cx[0]+r],[cx[1]-r,cx[1]+r],[cx[2]-r,cx[2]+r]]

def reimport(name):
    from importlib import reload
    reload(name)

def normalize(s):
    import string
    # 1. lowercase
    # 2. trim whitespaces
    for p in string.punctuation:
        s = s.replace(p, '')

    return s.lower().strip()

def fuzzymatch(s, answer_list, threshold=2):
    import string
    """
    Returns an answer among a list of answers that best matches given string.
    Maximum levenshtein distance is 2 by default.

    Refer to http://streamhacker.com/2011/10/31/fuzzy-string-matching-python/
    """
    for answer in answer_list:
        if normalize(answer) == normalize(s):
            return answer
    # no answer so far?

    edit_distance=[]
    for answer in answer_list:
        edit_distance.append(levenshtein(s, answer))

    min_ed = min(edit_distance)
    if min_ed <= threshold:
        return answer_list[edit_distance.index(min_ed)]
    else:
        print("No good matching. \
               The best one is {0} with the distance {1}".format(
               answer_list[edit_distance.index(min_ed)], mined))

def levenshtein(s1, s2):
    """
    Returns Levenshtein distance of two strings.

    levenshtein(s1, s2)

    Source: http://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance
    """
    if len(s1) < len(s2):
        return levenshtein(s2, s1)

    # len(s1) >= len(s2)
    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1 # j+1 instead of j since previous_row and current_row are one character longer
            deletions = current_row[j] + 1       # than s2
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row

    return previous_row[-1]



def dgyr2dnout(dt, nout_now):
    """
       returns nout closest to look back time at (nout) + dt.
       dt can be either positive or negative.

       example
       -------

       >>> nout_now = 180
       >>> nout_1gyr_earlier = dgyr2dnout(-1.0, nout_now)
       >>> nout_1gyr_earlier
           166

    """
    import general.defaults
    import numpy as np
    df = general.defaults.Default()
    df.load_time()

    lbt = df.times["lbt"]
    lbt_now = df.times["lbt"][df.times["nout"] == nout_now]
    #print(np.abs(lbt - lbt_now - dt))
    i_dt= np.argmin(np.abs(lbt - lbt_now + dt))

    return df.times["nout"][i_dt]


def get_last_snapshot(base='./'):
    """
    Some files/directories other than output_00xxx can be found under snapshots/
    In such case, listing everything under snapshots/ and sorting and selectnig 
    the last entry will fail (They may not even be a number.)

    Thus, search for output_ pattern and then search for the latest one.
    """

    fi = glob(base+"snapshots/output_?????/info*.txt")
    return int(np.sort(fi)[-1].split("info_")[1].split(".txt")[0])
    #return int(np.sort(os.listdir(base + "snapshots/"))[-1].split("_")[1])
