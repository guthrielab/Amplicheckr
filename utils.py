
#####Assorted helper functions#####

#translate a string into its complement
def complement(seq, bdict):
    compl = str()
    for i in str(seq):
        compl += bdict[i]
    return compl

#reverse a string
def reverse(forward):
    return forward[::-1]

#remove "-" from string
def removed(s):
    return (str(s).replace("-", ""))

#create the alignment visual with sticks based on 2 aligned seq
def alignvis (qalign, dalign, bnmatr):
    vislist = []
    for i, j in zip(qalign, dalign):
        if score(i,j, bnmatr)>0:
            vislist.append("|")
        else:
            vislist.append(" ")
    return "".join(vislist)

#return the score of two nt using a matrix
def score (kp, kt, matrix, convertd):
    return int(matrix[convertd[kp]][convertd[kt]])

def truncate (s):
    return s if len(s)<30 else s[:30]+"..."
