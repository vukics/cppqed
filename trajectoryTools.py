import numpy as np
import matplotlib.pyplot as plt
import json


# load trajectory
def loadTraj(fileDesc) :
    return np.loadtxt([l.replace('(',' ').replace(')',' ').replace(',',' ') for l in open(fileDesc,'r').readlines() if not l[0]=='#'])


# load jump trajectory
def jumpTraj(fileDesc) :
    return np.array(json.loads(open(fileDesc,'r').read().split("OUTRO:")[-1])['QJMC']['Jump trajectory'])


# photonClicks=(trajQJMC[trajQJMC[:,1]==1])[:,0]
# waitingTimes=photonClicks[1:]-photonClicks[:-1]


def concatenateByTime(trajs) :
    for i in range(1,len(trajs)) : trajs[i][:,0]+=trajs[i-1][-1,0]
    return np.concatenate(trajs)


def purgeRelaxation(traj,upper,i=6) :
    startingIndex=int(np.nonzero(traj[:,i]>upper)[0][0]/2)
    startingTime=traj[startingIndex,0]
    traj[startingIndex:,0]-=startingTime
    return traj[startingIndex:,...]


def relativeCorrelation(one,two) :
    return np.mean(one[:]*two[:])/(np.mean(one[:])*np.mean(two[:]))-1


def correlationInTime(traj) :
    return [relativeCorrelation(traj[i:],traj[:-i]) for i in range(100,100000,100)]


def upperPhotonNum(g) :
    return g**2/50.
