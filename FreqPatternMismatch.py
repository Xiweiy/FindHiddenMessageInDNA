#import sys
#lines = sys.stdin.read().splitlines()
#pattern = lines[0]
#dnaseq = lines[1]
#d = int(lines[2])

##Find out the most frequent pattern with mismatch and reverse complement seq

##Find out the min skew
def min_skew(dnaseq):
      GC = 0
      skew = [GC]
      for i in dnaseq:
            if i == 'G':
                  GC +=1
            elif i == 'C':
                  GC -=1
            skew.append(GC)

      import matplotlib.pyplot as plt
      plt.plot(skew)
      plt.savefig('skew.png')

      mincount = min(skew)
      for i in range(len(skew)):
            if skew[i] == mincount:
                  minposition = i
                  print i,

def count_Hamming_distance(seq1,seq2):
      count = 0
      for i,j in zip(seq1,seq2):
            if i !=j:
                  count = count +1
      return count


def neighbors(pattern, d):
      if not d:
            return pattern
      if len(pattern) ==1:
            return ['A', 'C', 'T', 'G']
      Neighborhood = []
      prefix = pattern[0]
      suffix = pattern[1:]
      suffixngbr = neighbors(suffix, d)
      for i in suffixngbr:
            if count_Hamming_distance(suffix, i) < d:
                  for j in ['A', 'C', 'T', 'G']:
                        Neighborhood.append(j+i)
            else:
                  Neighborhood.append(prefix+i)
      return Neighborhood


NucleoToNum = {'A':0, 'C':1, 'G':2, 'T':3}
def PatternToNumber(pattern):
    if len(pattern) ==1:
        return NucleoToNum[pattern]
    prefix = pattern[:len(pattern)-1]
    nucleo = pattern[-1]
    #print prefix, nucleo
    return 4*PatternToNumber(prefix) + NucleoToNum[nucleo]

NumToNucleo = {0:'A', 1:'C', 2:'G', 3:'T'}
def NumberToPattern(number, kmer):
    if kmer ==1:
        return NumToNucleo[number]
    quotient = number/4
    reminder = number%4
    return NumberToPattern(quotient, kmer-1) + NumToNucleo[reminder]

reverse_complementary={'A':'T', 'T':'A','G':'C','C':'G'}

def rev_comp(dnasequence):
      reverse_complement_sequence=""
      for i in dnasequence[::-1]:
          reverse_complement_sequence+=reverse_complementary[i]
      return reverse_complement_sequence

def FreqWordMismatch_sort(dnaseq, k, d):
      neighborhood = []
      nseq = len(dnaseq)
      for i in range(nseq-k+1):
            pattern = dnaseq[i:i+k]
            neighborhood += neighbors(pattern, d)
      index = []
      for i in neighborhood:
            index.append(PatternToNumber(i))
      n = len(index)
      count = [1]*n
      index.sort()

      maxcount = 0
      FreqPattern = []
      for i in range(1,n):
            if index[i] == index[i-1]:
                  count[i] = count[i-1] +1
      maxcount = max(count)
      for i in range(n):
            if count[i] == maxcount:
                  FreqPattern.append(NumberToPattern(index[i],k))
      return FreqPattern

#FreqPattern = FreqWordMismatch_sort(dnaseq, k,d)
#for i in FreqPattern:
#      print i,
#print

def FreqWordMismatchRevComp(dnaseq, k ,d):
      neighborhood = []
      nseq = len(dnaseq)
      for i in range(nseq-k+1):
            pattern = dnaseq[i:i+k]
            neighborhood += neighbors(pattern, d)
            reverse_pattern = rev_comp(pattern)
            neighborhood += neighbors(reverse_pattern, d)
      index = []
      for i in neighborhood:
            index.append(PatternToNumber(i))
      n = len(index)
      count = [1]*n
      index.sort()

      maxcount = 0
      FreqPattern = []
      for i in range(1,n):
            if index[i] == index[i-1]:
                  count[i] = count[i-1] +1
      maxcount = max(count)
      for i in range(n):
            if count[i] == maxcount:
                  FreqPattern.append(NumberToPattern(index[i],k))
      return FreqPattern


FreqPattern = FreqWordMismatchRevComp(dnaseq, k,d)
for i in FreqPattern:
      print i,

