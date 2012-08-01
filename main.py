#! /usr/bin/env python

import math
import numpy as np
import re
import sys

from optparse import OptionParser

class RotationMatrix(np.ndarray):

    def __new__(self, angle):
        return np.array([[round(math.cos(angle)), round(-math.sin(angle))],
            [round(math.sin(angle)), round(math.cos(angle))]])


class Rotation(object):

    def __init__(self, point):
        self.point = point
    
    def rotate(self):
        pass


class R90(Rotation):

    matrix = RotationMatrix(0.5*math.pi)


class R180(Rotation):

    matrix = RotationMatrix(math.pi)


class R270(Rotation):

    matrix = RotationMatrix(1.5*math.pi)


class Aminoacid(str):

    HYDROPHOBIC = 'H'
    POLAR = 'P'
    
    def isHydrophobic(self):
        return self is self.HYDROPHOBIC


class Chain(str):
    
    def __getitem__(self, index):
        return Aminoacid(str.__getitem__(self, index))
        
    def __new__(self, value): 
        if not re.match('[HP]+$', value):
            raise ValueError("Sequence should consist of 'H' and 'P' characters")
        return str.__new__(self, value)


class Microstate(object):

    chain = ''
    
    def __init__(self, coords, energy=10000):
        self.coords = coords
        self.energy = energy
    
    def setEnergy(self):
        pass

    def transform(self):
        pass


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
    
    Microstate.chain = Chain(options.sequence)


if __name__ == '__main__':
    main()
