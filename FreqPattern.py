##The Efficient algm to find freq pattern in a genome using the FreqArray Approach


##Get standard input
"""import sys 
lines = sys.stdin.read().splitlines()
dnasequence = lines[0]
k, L, t = lines[1].split(" ")
k=int(k)
L=int(L)
t=int(t)"""

#infile = open('E-coli.txt','r')
#inlist = infile.readlines()
#dnasequence = inlist[0].rstrip('\n')
#k, L, t = inlist[1].rstrip('\n').split(" ")

dnasequence  ='ACGTACGT'
k=1
L=5
t=2
##Initiate a FreqArray to store the occurance of a pattern in a given window
##a Clump array to store whether a pattern has ever been frequent in any window
##an empty freq pattern
FreqArray = [0]*4**k
clump = [0]*4**k
FreqPattern = []

##Convert pattern to position(number)
NucleoToNum = {'A':0, 'C':1, 'G':2, 'T':3}
def PatternToNumber(pattern):
    if len(pattern) ==1:
        return NucleoToNum[pattern]
    prefix = pattern[:len(pattern)-1]
    nucleo = pattern[-1]
    #print prefix, nucleo
    return 4*PatternToNumber(prefix) + NucleoToNum[nucleo]

##Convert position(number) to pattern
NumToNucleo = {0:'A', 1:'C', 2:'G', 3:'T'}
def NumberToPattern(number, kmer):
    if kmer ==1:
        return NumToNucleo[number]
    quotient = number/4
    reminder = number%4
    return NumberToPattern(quotient, kmer-1) + NumToNucleo[reminder]

##Compute the # of occurance of a particular pattern, return the frequency array
def CountFreq(string, kmer):
    for i in range(len(string)-kmer+1):
        pattern = string[i: i+kmer]
        position = PatternToNumber(pattern)
        #print pattern, position
        FreqArray[position] += 1
        #print FreqArray
    maxcount = max(FreqArray)
    for i in range(4**k):
        if FreqArray[i] >= t:
            clump[i] = 1


CountFreq(dnasequence[0:0+L], k)
for i in range(1,len(dnasequence)-L+1):
    firstpattern = dnasequence[i-1:i+k-1]
    FreqArray[PatternToNumber(firstpattern)] -= 1
    lastpattern = dnasequence[i+L-k: L+i]
    lastpos = PatternToNumber(lastpattern)
    FreqArray[lastpos] +=1
    if FreqArray[lastpos] >= t:
        clump[lastpos] = 1
for i in range(4**k):
    if clump[i] ==1:
        FreqPattern.append(NumberToPattern(i, k))
        print NumberToPattern(i, k),

#print len(FreqPattern)
