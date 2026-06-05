#####Global dictionaries for scoring, conversion, and complement (ran once only)#####
import numpy as np

#score matrix used in blastn/initial alignment
def bnmatrix():
    '''rawmat = """
    A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N
A   2  -3  -3  -3  -3  -1  -1  -3  -3  -1  -3  -1  -1  -1  -2
T  -3   2  -3  -3  -3  -1  -3  -1  -1  -3  -1  -3  -1  -1  -2
G  -3  -3   2  -3  -1  -3  -1  -3  -1  -3  -1  -1  -3  -1  -2
C  -3  -3  -3   2  -1  -3  -3  -1  -3  -1  -1  -1  -1  -3  -2
S  -3  -3  -1  -1  -1  -3  -1  -1  -1  -1  -1  -1  -1  -1  -2
W  -1  -1  -3  -3  -3  -1  -1  -1  -1  -1  -1  -1  -1  -1  -2
R  -1  -3  -1  -3  -1  -1  -1  -3  -1  -1  -1  -1  -1  -1  -2
Y  -3  -1  -3  -1  -1  -1  -3  -1  -1  -1  -1  -1  -1  -1  -2
K  -3  -1  -1  -3  -1  -1  -1  -1  -1  -3  -1  -1  -1  -1  -2
M  -1  -3  -3  -1  -1  -1  -1  -1  -3  -1  -1  -1  -1  -1  -2
B  -3  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -2
V  -1  -3  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -2
H  -1  -1  -3  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -2
D  -1  -1  -1  -3  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -2
N  -2  -2  -2  -2  -2  -2  -2  -2  -2  -2  -2  -2  -2  -2  -2""" '''
    rawmat= """
    A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N
A   2  -3  -3  -3  -3   2   2  -3  -3   2  -3   2   2   2   2
T  -3   2  -3  -3  -3   2  -3   2   2  -3   2  -3   2   2   2
G  -3  -3   2  -3   2  -3   2  -3   2  -3   2   2  -3   2   2
C  -3  -3  -3   2   2  -3  -3   2  -3   2   2   2   2  -3   2
S  -3  -3   2   2   2  -3   2   2   2   2   2   2   2   2   2
W   2   2  -3  -3  -3   2   2   2   2   2   2   2   2   2   2
R   2  -3   2  -3   2   2   2  -3   2   2   2   2   2   2   2
Y  -3   2  -3   2   2   2  -3   2   2   2   2   2   2   2   2
K  -3   2   2  -3   2   2   2   2   2  -3   2   2   2   2   2
M   2  -3  -3   2   2   2   2   2  -3   2   2   2   2   2   2
B  -3   2   2   2   2   2   2   2   2   2   2   2   2   2   2
V   2  -3   2   2   2   2   2   2   2   2   2   2   2   2   2
H   2   2  -3   2   2   2   2   2   2   2   2   2   2   2   2
D   2   2   2  -3   2   2   2   2   2   2   2   2   2   2   2
N   2   2   2   2   2   2   2   2   2   2   2   2   2   2   2"""
    
    removealpha = rawmat.strip().split()
    abcd = [i for i in removealpha if not i.isalpha()]
    
    blastn = np.array(abcd).reshape(15, 15)
    
    #above translation in np array (0-14 by 0-14)
    return blastn

#score matrix used for score the mismatches between primer and its complement genome
#note pur/pur is worse than py/py is worse than pur/py
def mmmatrix():
    rawmat= """
    A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N
A  -5  +2  -5  -3  -3  +2  -5  +2  +2  -3  +2  +2  +2  -3  +2
T  +2  -4  -3  -4  -3  +2  +2  -4  -3  +2  -3  +2  +2  +2  +2
G  -5  -3  -5  +2  +2  -3  -5  +2  -3  +2  +2  -3  +2  +2  +2
C  -3  -4  +2  -4  +2  -3  +2  -4  +2  -3  +2  +2  -3  +2  +2
S  -3  -3  +2  +2  +2  -3  +2  +2  +2  +2  +2  +2  +2  +2  +2
W  +2  +2  -3  -3  -3  +2  +2  +2  +2  +2  +2  +2  +2  +2  +2
R  -5  +2  -5  +2  +2  +2  -5  +2  +2  +2  +2  +2  +2  +2  +2
Y  +2  -4  +2  -4  +2  +2  +2  -4  +2  +2  +2  +2  +2  +2  +2
K  +2  -3  -3  +2  +2  +2  +2  +2  -3  +2  +2  +2  +2  +2  +2
M  -3  +2  +2  -3  +2  +2  +2  +2  +2  -3  +2  +2  +2  +2  +2
B  +2  -3  +2  +2  +2  +2  +2  +2  +2  +2  +2  +2  +2  +2  +2
V  -3  +2  +2  +2  +2  +2  +2  +2  +2  +2  +2  +2  +2  +2  +2
H  +2  +2  +2  -3  +2  +2  +2  +2  +2  +2  +2  +2  +2  +2  +2
D  +2  +2  -3  +2  +2  +2  +2  +2  +2  +2  +2  +2  +2  +2  +2
N  +2  +2  +2  +2  +2  +2  +2  +2  +2  +2  +2  +2  +2  +2  +2
"""
    removealpha = rawmat.strip().split()
    abcd = [i for i in removealpha if not i.isalpha()]
    
    mmatrix = np.array(abcd).reshape(15, 15)
    
    #above translation in np array (0-14 by 0-14)
    return mmatrix

#conversion of nt value to 0-14 for blastn matrix use (ex. G + T mm -> blastn[2][1] = -3)
def convertdict ():
    rawline = "A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N"
    ntbase = rawline.split()
    conversiondict = dict((i, ntbase.index(i)) for i in ntbase)
    
    #dict of A: 0, T: 1, etc
    return conversiondict

#dict of A map to T, T map to A for complement use
def basedict ():
    basedict = {"A": "T", "C": "G", "R": "Y", "W": "S", "M": "K", "V": "B", "H": "D"}
    basedict.update(dict((basedict[k], k) for k in basedict))
    basedict["N"] = "N"
    #dict of A: T, T:A, C:G, G:C, etc.
    return basedict
