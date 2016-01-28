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

def MOTIFENUMERATION(DNA, K, D):
      pattern = {}
      for each_string in DNA:
            for i in range(len(each_string)-K+1):
                  kmer = each_string[i:i+K]
                  kmer_ngbr = neighbors(kmer, D)
                  if not D:
                        kmer_ngbr = [kmer_ngbr]
                  for each_ngbr in kmer_ngbr:
                        counter = 0
                        for each_dna in DNA:
                              count = 0
                              for j in range(len(each_dna)-K+1):
                                    kmer2 = each_dna[j:j+K]
                                    if count_Hamming_distance(each_ngbr, kmer2) <=D:
                                          count += 1
                              if not count:
                                    break
                              counter +=1
                        if counter == len(DNA):
                              pattern[each_ngbr] = 0
      return pattern.keys()
                        
def DISTANCEPATTERNSTRING(PATTERN,DNA):
      k = len(PATTERN)
      distance = 0
      for each_string in DNA:
            Hamming = 1000
            for i in range(len(each_string)-k+1):
                  pattern = each_string[i:i+k]
                  d = count_Hamming_distance(pattern, PATTERN)
                  if d < Hamming:
                        Hamming = d
            distance += Hamming
      return distance

#import sys
#lines = sys.stdin.read().splitlines()
#pattern = lines[0].rstrip('\n')
#DNA = lines[1].rstrip('\n')
#DNA = DNA.split(' ')
#print DISTANCEPATTERNSTRING(pattern, DNA)

NumToNucleo = {0:'A', 1:'C', 2:'G', 3:'T'}
def NumberToPattern(number, kmer):
    if kmer ==1:
        return NumToNucleo[number]
    quotient = number/4
    reminder = number%4
    return NumberToPattern(quotient, kmer-1) + NumToNucleo[reminder]

def MEDIANSTRING(K, DNA):
      distance = 1000000
      for i in range(4**K):
            pattern = NumberToPattern(i, K)
            d = DISTANCEPATTERNSTRING(pattern, DNA)
            if d < distance:
                  distance = d
                  median = pattern
      return median

#import sys
#lines = sys.stdin.read().splitlines()
#K = int(lines[0].rstrip('\n'))
#DNA = [i.rstrip('\n') for i in lines[1:]]
#print MEDIANSTRING(K, DNA)

from numpy import array

rowindex = {'A':0, 'C':1, 'G':2, 'T':3}
def MOSTPROBABLE(DNA, K, profile):
      maxprob = 0
      PATTERN = DNA[:K]
      for i in range(len(DNA)-K+1):
            pattern = DNA[i:i+K]
            prob = 1
            for i in range(K):
                  prob = prob * profile[rowindex[pattern[i]],i]
            if prob > maxprob:
                  maxprob = prob
                  PATTERN = pattern
      return PATTERN

#import sys
#lines = sys.stdin.read().splitlines()
#DNA = lines[0].rstrip('\n')
#K = int(lines[1].rstrip('\n'))
#profile = ''
#for i in lines[2:]:
#      profile = profile + i.rstrip('\n') + ' '
#profile = array([float(i) for i in profile.rstrip(' ').split(' ')]).reshape((4, K))


from pandas import DataFrame
from copy import deepcopy

def UPDATEPROFILE(OLDMATRIX, NEW_MOTIF):
      for i in range(len(NEW_MOTIF)):
            OLDMATRIX[rowindex[NEW_MOTIF[i]],i] +=1
      return OLDMATRIX
      
def FINDCONCENSUS(profile):
      concensus = ''
      profile = DataFrame(profile)
      for i in range(K):
            maxprob = max(profile[i])
            row = profile.index[profile[i]==maxprob][0]
            concensus += revrowindex[row]
      return concensus

revrowindex = {0:'A', 1:'C', 2:'G',3:'T'}
def GREEDYMOTIF(DNA,K,T):
      BESTMOTIF = []
      countmatrix = array([[1.0]*K]*4)
      for i in DNA:
            BESTMOTIF.append(i[:K])
            countmatrix = UPDATEPROFILE(countmatrix, BESTMOTIF[-1])
      profile = countmatrix/len(BESTMOTIF)
      concensus = FINDCONCENSUS(profile)
      bestscore = DISTANCEPATTERNSTRING(concensus, BESTMOTIF)
      for i in range(len(DNA[0])-K+1):
            MOTIF = [DNA[0][i:i+K]]
            #print MOTIF[0]
            countmatrix = UPDATEPROFILE(array([[1.0]*K]*4), MOTIF[0])
            #print countmatrix
            profile = countmatrix
            for each_dna in DNA[1:]:
                  MOTIF.append(MOSTPROBABLE(each_dna, K, profile))
                  countmatrix = UPDATEPROFILE(countmatrix, MOTIF[-1])
                  profile = countmatrix/len(MOTIF)
            concensus = FINDCONCENSUS(profile)
            curscore = DISTANCEPATTERNSTRING(concensus, MOTIF)
            if curscore < bestscore:
                  BESTMOTIF = deepcopy(MOTIF)
                  bestscore = curscore
      return BESTMOTIF

#import sys
#lines = sys.stdin.read().splitlines()
#K, T= lines[0].rstrip('\n').split(' ')
#K = int(K)
#T = int(T)
#DNA = [i.rstrip('\n') for i in lines[1:]]

#MOTIFS = GREEDYMOTIF(DNA, K,T)
#for i in MOTIFS:
#      print i
