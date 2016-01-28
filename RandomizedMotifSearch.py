#import sys
#lines = sys.stdin.read().splitlines()
#K, T = lines[0].rstrip('\n').split(' ')
#K = int(K)
#T = int(T)
#DNA = [i.rstrip('\n') for i in lines[1:]]

from pandas import DataFrame
from numpy import random, array
from copy import deepcopy

random.seed(1)
def count_Hamming_distance(seq1,seq2):
      count = 0
      for i,j in zip(seq1,seq2):
            if i !=j:
                  count = count +1
      return count

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
def RANDOMIZEDMOTIF(DNA,K,T):
      MOTIF = []
      randn = random.random_integers(0, len(DNA[0]), T)
      countmatrix = array([[1.0]*K]*4)
      for i in range(T):
            index = randn[i]
            MOTIF.append(DNA[i][index:index+K])
            countmatrix = UPDATEPROFILE(countmatrix, MOTIF[-1])

      curscore = 9999
      bestscore = 10000
      while curscore < bestscore:
            bestscore = curscore

            profile = countmatrix/T
            concensus = FINDCONCENSUS(profile)
            curscore = DISTANCEPATTERNSTRING(concensus, MOTIF)
            
            MOTIF = []
            countmatrix = array([[1.0]*K]*4)
            for each_dna in DNA:
                  MOTIF.append(MOSTPROBABLE(each_dna, K, profile))
                  countmatrix = UPDATEPROFILE(countmatrix, MOTIF[-1])
      return MOTIF, bestscore

'''
bestscore = 1000
for i in range(1000):
      motif, best = RANDOMIZEDMOTIF(DNA, K, T)      
      if best < bestscore:
            MOTIF = motif
            bestscore = best
            #print bestscore, MOTIF
for i in MOTIF:
      print i
'''


K = 8
T = 5
N = 100
DNA = ['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA',
'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
'''
import sys
lines = sys.stdin.read().splitlines()
K, T, N = lines[0].rstrip('\n').split(' ')
K = int(K)
T = int(T)
N = int(N)
DNA = [i.rstrip('\n') for i in lines[1:]]
'''

def UPDATEPROFILE2(OLDMATRIX, NEW_MOTIF):
      for i in range(len(NEW_MOTIF)):
            OLDMATRIX[rowindex[NEW_MOTIF[i]],i] -=1
      return OLDMATRIX

def FINDPROBS(DNA, K, profile):
      probs = []
      for i in range(len(DNA)-K+1):
            pattern = DNA[i:i+K]
            prob = 1
            for i in range(K):
                  prob = prob * profile[rowindex[pattern[i]],i]
            probs.append(prob)
      return array(probs)

      
def GIBBSSAMPLER(DNA, K, T, N):
      MOTIF = []
      randn = random.random_integers(0, len(DNA[0]), T)
      countmatrix = array([[1.0]*K]*4)
      for i in range(T):
            index = randn[i]
            MOTIF.append(DNA[i][index:index+K])
            countmatrix = UPDATEPROFILE(countmatrix, MOTIF[-1])
      
      bestscore = 1000
      BESTGROUP = []
      for i in range(N):
            randdna = random.random_integers(0, T-1)
            countmatrix = UPDATEPROFILE2(countmatrix, MOTIF[randdna])  
            profile = countmatrix/(T-1)
            probs = FINDPROBS(DNA[randdna], K, profile)
            probs = probs/sum(probs)
            randi = random.choice(range(len(probs)), p = probs)
            motifi = DNA[randdna][randi:randi+K]
            MOTIF[randdna] = motifi

            countmatrix = UPDATEPROFILE(countmatrix, motifi)
            profile = countmatrix/T
            concensus = FINDCONCENSUS(profile)
            curscore = DISTANCEPATTERNSTRING(concensus, MOTIF)
            if curscore < bestscore:
                  BESTGROUP = deepcopy(MOTIF)
                  bestscore = curscore
      return BESTGROUP, bestscore
##Note: since list is copied by reference, then BESTGROUP will change if MOTIF change
##     have to make a deep copy

bestscore = 1000
for i in range(20):
      print i
      motif, score =  GIBBSSAMPLER(DNA, K, T, N)
      print motif, score
      if score < bestscore:
            MOTIF = motif
            bestscore = score
for i in MOTIF:
      print i

