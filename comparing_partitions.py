import sys
import math
import pandas 
import numpy as np

def getContTable(ar1, ar2):
    cont = {}
    for i in xrange(0, len(ar1)):
        keyAr1 = ar1[i]
        keyAr2 = ar2[i]
        if keyAr1 in cont:
            if keyAr2 in cont[keyAr1]:
                cont[keyAr1][keyAr2] += 1
            else:
                cont[keyAr1][keyAr2] = 1
        else:
            cont[keyAr1] = {keyAr2: 1}
    return cont


def getContTableTotals(cont, ar1, ar2):
    sumRow = {}
    sumCol = {}
    h1 = set(ar1)
    h2 = set(ar2)
    for x in h2:
        sumRow[x] = 0
        for y in h1:
            if y in cont:
                if x in cont[y]:
                    val = cont[y][x]
                    sumRow[x] += val
                    if y in sumCol:
                        sumCol[y] += val
                    else:
                        sumCol[y] = val
    total = 0
    for x in h1:
        total += sumCol[x]
    return (sumRow, sumCol, total)


def getMismatchMatrix(cont, ar1, ar2):
    totals = getContTableTotals(cont, ar1, ar2)
    # print totals
    n = totals[2]
    h1 = set(ar1)
    h2 = set(ar2)
    a = 0
    for x in h1:
        for y in h2:
            if x in cont:
                if y in cont[x]:
                    val = cont[x][y]
                    a += (val * (val - 1)) / 2
    a1 = 0
    sumCol = totals[1]
    for x in sumCol:
        val = sumCol[x]
        a1 += (val * (val - 1)) / 2

    b = a1 - a

    a2 = 0
    sumRow = totals[0]
    for x in sumRow:
        val = sumRow[x]
        a2 += (val * (val - 1)) / 2
    c = a2 - a
    d = float((n * (n - 1)) / 2) - a1 - c
    return (a, b, c, d, n)


def getRand(a, b, c, d):
    rand = (a + d) / float(a + b + c + d)
    return rand


def getWallace(a, b, c):
    w1 = float(0)
    w2 = float(0)
    if (a + b) > 0:
        w1 = a / float(a + b)
    if (a + c) > 0:
        w2 = a / float(a + c)
    return (w1, w2)


def getSimpsons(ar):
    n = len(ar)
    d = {}
    for x in ar:
        if x in d:
            d[x] += 1
        else:
            d[x] = 1

    sumTotal = 0
    sumFc2 = 0
    sumFc3 = 0
    for x in d:
        val = d[x]
        sumTotal += val * (val - 1)
        sumFc2 += (val / float(n)) ** 2.0
        sumFc3 += (val / float(n)) ** 3.0

    sid = 1.0
    if n * (n - 1) > 0:
        sid = 1.0 - (sumTotal / float(n * (n - 1)))

    sqSumFc2 = sumFc2 ** 2.0
    s2 = (4.0 / float(n)) * float(sumFc3 - sqSumFc2)

    sidLow = sid - 2 * math.sqrt(s2)
    sidHigh = sid + 2 * math.sqrt(s2)

    return (sid, sidLow, sidHigh, len(d))


def getAdjustedWallace(cont, ar1, ar2):
    sumCol, sumRow, totals = getContTableTotals(cont, ar1, ar2)
    rSumW1 = 0
    rSumW2 = 0
    csumFc2 = {}
    csumFc3 = {}
    for i in cont:
        rSum = sumRow[i]
        rsumFc2 = 0
        rsumFc3 = 0
        for j in cont[i]:
            val = cont[i][j]
            cSum = sumCol[j]
            rsumFc2 += (val / float(rSum)) ** 2.0
            rsumFc3 += (val / float(rSum)) ** 3.0
            if j in csumFc2:
                csumFc2[j] += (val / float(cSum)) ** 2.0
            else:
                csumFc2[j] = (val / float(cSum)) ** 2.0
            if j in csumFc3:
                csumFc3[j] += (val / float(cSum)) ** 3.0
            else:
                csumFc3[j] = (val / float(cSum)) ** 3.0
        rsqrsumFc2 = rsumFc2 ** 2.0
        rvarSID = 0.0
        if rSum > 1:
            rvarSID = float(4.0 * rSum * (rSum - 1.0) * (rSum - 2.0) * rsumFc3 + 2.0 * rSum * (rSum - 1.0) * rsumFc2 - 2.0 * rSum * (rSum - 1.0) * (2.0 * rSum - 3.0) * rsqrsumFc2) / float((rSum * (rSum - 1.0)) ** 2.0)
        rSumW1 += float((rSum * (rSum - 1.0)) ** 2.0) * rvarSID
        rSumW2 += rSum * (rSum - 1.0)

    csumw1 = 0
    csumW2 = 0

    for j in sumCol:
        cSum = sumCol[j]
        fc2 = csumFc2[j]
        fc3 = csumFc3[j]
        csqrsumFc2 = fc2 ** 2.0
        cvarSID = 0.0
        if cSum > 1:
            # (4*sumcol[j]*(sumcol[j]-1)*(sumcol[j]-2)*csumFc3[j]+2*sumcol[j]*(sumcol[j]-1)*csumFc2[j]-2*sumcol[j]*(sumcol[j]-1)*(2*sumcol[j]-3)*csqrsumFc2[j])/sqr(sumcol[j]*(sumcol[j]-1));
            cvarSID = float(4.0 * cSum * (cSum - 1.0) * (cSum - 2.0) * fc3 + 2.0 * cSum * (cSum - 1.0) * fc2 - 2.0 * cSum * (cSum - 1.0) * (2.0 * cSum - 3.0) * csqrsumFc2) / float((cSum * (cSum - 1.0)) ** 2.0)
            #print cvarSID
        csumw1 += (cSum * (cSum - 1.0) * cvarSID) ** 2.0
        csumW2 += cSum * (cSum - 1.0)
    #print "csumW2=%f" % csumW2
    varW1 = 0.0
    varW2 = 0.0
    if rSumW2 > 0:
        varW1 = float(rSumW1 / float(float(rSumW2) ** 2.0))
    if csumW2 > 0:
        varW2 = float(csumw1 / float(float(csumW2) ** 2.0))
    #print "varW1=%f" % varW1
    #print "varW2=%f" % varW2
    a, b, c, d, n = getMismatchMatrix(cont, ar1, ar2)
    w1, w2 = getWallace(a, b, c)
    sid1 = getSimpsons(ar1)
    sid2 = getSimpsons(ar2)
    wi1 = 1 - sid1[0]
    wi2 = 1 - sid2[0]
    aw1 = float(w1 - wi2) / float(1 - wi2)
    aw2 = float(w2 - wi1) / float(1 - wi1)
    
    aw1CI = 2.0 * (1.0 / float(1.0 - wi2)) * math.sqrt(varW1)
    aw1Low = aw1 - aw1CI
    aw1High = aw1 + aw1CI
    
    aw2CI = 2.0 * (1.0 / float(1.0 - wi1)) * math.sqrt(varW2)
    #print "wi1=%f" % wi1
    #print "aw2CI=%f" % aw2CI
    aw2Low = aw2 - aw2CI
    aw2High = aw2 + aw2CI
    return (aw1, aw1Low, aw1High, aw2, aw2Low, aw2High)


print len(sys.argv)  
if len(sys.argv) == 2:
    filename  = sys.argv[1]
    data = pandas.read_table(filename)
    
    ncols = len(data.columns) 

    for col1 in data.columns:
        #print i
        for col2 in data.columns:
            #print j
            if col1 != col2:
                print "===========",col1,"=",col2,"============"
                #print data[col2]
                icol = data[col1].values.tolist()
                jcol = data[col2].values.tolist()
                # print y
                #print "================="
                # print icol
                # print "================="
                # print jcol
                # print "================="

                cont = getContTable(icol, jcol)

                #print cont
                #sys.exit()
                abcdn = getMismatchMatrix(cont, icol, jcol)
                #print "Rand\t%f" % getRand(abcdn[0], abcdn[1], abcdn[2], abcdn[3])
                wallace = getWallace(abcdn[0], abcdn[1], abcdn[2])
                print "Wallace 1vs2\t%f" % wallace[0]
                print "Wallace 2vs1\t%f" % wallace[1]
                sid1 = getSimpsons(icol)
                sid2 = getSimpsons(jcol)
                #print sid1
                print "Simpsons 1\t%f\t%f\t%f\t%d" % sid1
                print "Simpsons 2\t%f\t%f\t%f\t%d" % sid2
                #sys.exit()
                AdjustedWallace = getAdjustedWallace(cont, icol, jcol)
                print "Adjusted Wallace 1vs2\t%f\t%f\t%f" % AdjustedWallace[0:3]
                print "Adjusted Wallace 2vs1\t%f\t%f\t%f" % AdjustedWallace[3:6] # STILL BROKEN #FIXIT

#else:
#    print """
# Please specify target file for comparisons

# Output:

# Simpson's index of diversity
# Wallace coefficient 1vs2
# Wallace coefficient 2vs1
# Adjusted Wallace coefficient 1vs2
# Adjusted wallace coefficient 2vs1


    """


# if len(sys.argv) == 3:
#     array1 = sys.argv[1]
#     array2 = sys.argv[2]
#     x = array1.split(",")
#     # print x
#     y = array2.split(",")
#     # print y
#     cont = getContTable(x, y)
#     abcdn = getMismatchMatrix(cont, x, y)
#     print "Rand\t%f" % getRand(abcdn[0], abcdn[1], abcdn[2], abcdn[3])
#     wallace = getWallace(abcdn[0], abcdn[1], abcdn[2])
#     print "Wallace 1vs2\t%f" % wallace[0]
#     print "Wallace 2vs1\t%f" % wallace[1]
#     sid1 = getSimpsons(x)
#     sid2 = getSimpsons(y)
#     print "Simpsons 1\t%f\t%f\t%f\t%d" % sid1
#     print "Simpsons 2\t%f\t%f\t%f\t%d" % sid2
#     print getAdjustedWallace(cont, x, y)
# else:
#     print """
# Please specify 2 lists of clusters delimited by commas ','
# Output:
# Rand coefficient
# Wallace coefficient 1vs2
# Wallace coefficient 2vs1
# Simpson's index of diversity
#     """
