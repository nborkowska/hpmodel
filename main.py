#! /usr/bin/env python

#
# python 2.7
#

import math
import numpy as np
import random
import re
import sys

from optparse import OptionParser

class RotationMatrix(np.ndarray):

    def __new__(self, angle):
        return np.array([[round(math.cos(angle)), round(-math.sin(angle))],
            [round(math.sin(angle)), round(math.cos(angle))]])


class Rotation(object):

    matrices = [RotationMatrix(angle) \
            for angle in [0.5*math.pi, math.pi, 1.5*math.pi]]
    
    def __init__(self, point):
        self.point = point
    
    def rotate(self, inputDict):           #TODO improve
        matrix = random.choice(self.matrices)
        index = inputDict.get(self.point).id
        newCoords = {}
        for key in inputDict:
            aminoacid = inputDict.get(key)
            if aminoacid.id > index:
                key = tuple(np.add(self.point,\
                        np.dot(matrix,np.subtract(key,self.point))))
            if not newCoords.get(key, 0):
                newCoords[key] = aminoacid
            else:
                return {}
        return newCoords


class Aminoacid(object):

    HYDROPHOBIC = 'H'
    POLAR = 'P'

    def __init__(self, id, symbol):
        self.id = id
        self.symbol = symbol


class Chain(list):

    def __new__(self, sequence): 
        if not re.match('[HP]+$', sequence):
            raise ValueError("Sequence should consist of 'H' and 'P' characters")
        return [Aminoacid(index,char) for index, char in enumerate(sequence)]


class Microstate(object):
    
    def __init__(self, coords):
        self.coords = coords
   
    def calculateEnergy(self):
        pass 

    def transform(self):
        """ returns next microstate """
        newCoords = {}
        while not newCoords:
            rotationPoint = random.choice(self.coords.keys())
            newCoords = Rotation(rotationPoint).rotate(self.coords)
        return Microstate(newCoords)
    
    @staticmethod
    def getInitialCoords(chain):
        return {(len(chain)/2,v):k for v,k in enumerate(chain)}


class Metropolis(object):
    
    def metropolis():
        pass


class SimulatedAnnealing(object):
    
    def __init__(self, **kwargs):
        kwargs.update({'chain': Chain(kwargs.pop('sequence'))})
        self.__dict__.update(kwargs)


def main():
    
    parser = OptionParser(usage="usage: %prog [options]", version="%prog 1.0")

    parser.add_option('-s', '--sequence', dest='sequence', \
            default='', help='protein sequence, ex. "HPHPHPHPPHPH"')
    parser.add_option('-i', '--maximum', type='float', dest='tMax', \
            default=1.0, help='initial temperature')
    parser.add_option('-f', '--minimum', type='float', dest='tMin', \
            default=0.15, help='final temperature')
    parser.add_option('-d', '--delta', type='float', dest='delta', \
            default=0.05, help='temperature decrement factor')
    parser.add_option('-k', '--transitions', type='int', dest='noTransitions',\
            default=10000, help='the length of the k-th Markov chain')

    (options, args) = parser.parse_args()
    
    simulation = SimulatedAnnealing(**options.__dict__)
    b=Metropolis()
    c = Microstate.getInitialCoords(simulation.chain)
    #print c, type(c)
    d = Microstate(c)
    d.transform()
    print d.coords

if __name__ == '__main__':
    main()
