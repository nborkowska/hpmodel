#! /usr/bin/env python

#
# python 2.7
#

import math
import numpy as np
import random
import re

from matplotlib import pyplot as plt
from matplotlib import cm
from optparse import OptionParser


def main():
    
    parser = OptionParser(usage="usage: %prog [options]", version="%prog 1.0")

    parser.add_option('-s', '--sequence', dest='sequence', \
            default='', help='protein sequence, ex. "PHPPHPPHHPPHHPPHPPHP"')
    parser.add_option('-i', '--maximum', type='float', dest='tMax', \
            default=1.0, help='initial temperature')
    parser.add_option('-f', '--minimum', type='float', dest='tMin', \
            default=0.15, help='final temperature')
    parser.add_option('-d', '--delta', type='float', dest='delta', \
            default=0.05, help='temperature decrement factor')
    parser.add_option('-k', '--transitions', type='int', dest='noTransitions',\
            default=10, help='the length of the k-th Markov chain')
    parser.add_option('-r', '--replicas', type='int', dest='replicas',\
            default=5, help='no of replicas')
    parser.add_option('-l', '--length', type='int', dest='length',\
            default=100, help='simulation length')
    parser.add_option('-e', '--frequency', type='int', dest='frequency',\
            default=10, help='frequency of replica exchange')
    parser.add_option('-t', '--type', type='choice', dest='type',\
            choices=['sa', 're'], default='re', \
            help='type of simulation: simulated annealing or replica exchange')

    #TODO check if combination of args is appropriate
    (options, args) = parser.parse_args()
    
    if options.type == 're':
        simulation = ReplicaExchange(**options.__dict__)
    else:
        simulation = SimulatedAnnealing(**options.__dict__)
    
    simulation.start()

def getId():
    x = 1
    while True:
        yield x
        x +=1


class RotationMatrix(np.ndarray):

    def __new__(self, angle):
        return np.array([[round(math.cos(angle)), round(-math.sin(angle))],
            [round(math.sin(angle)), round(math.cos(angle))]])


class Rotation(object):

    matrices = [RotationMatrix(angle) \
            for angle in [0.5*math.pi, math.pi, 1.5*math.pi]]
    
    def __init__(self, point):
        self.point = point
    
    def rotate(self, inputDict):           #TODO
        matrix = random.choice(self.matrices)
        index = inputDict.get(self.point).Id
        newCoords = {}
        for key in inputDict:
            aminoacid = inputDict.get(key)
            if aminoacid.Id > index:
                key = tuple(np.add(self.point,\
                        np.dot(matrix,np.subtract(key,self.point))))
            if not newCoords.get(key, None):
                newCoords[key] = aminoacid
            else:
                return {}
        return newCoords


class Aminoacid(object):

    HYDROPHOBIC = 'H'
    POLAR = 'P'

    def __init__(self, Id, symbol):
        self.Id = Id
        self.symbol = symbol


class Chain(list):

    def __new__(self, sequence): 
        if not re.match('[HP]+$', sequence):
            raise ValueError("Sequence should consist of 'H' and 'P' characters")
        return [Aminoacid(index,char) for index, char in enumerate(sequence)]


class Microstate(object):
    
    directions = [(1,0),(0,1)]
    ids = getId()

    def __init__(self, coords, energy=0.0, noContacts=0):
        self.Id = self.ids.next() 
        self.coords = coords
        self.energy = energy
        self.contacts = noContacts
   
    def setContacts(self, noContacts):
        self.contacts = noContacts

    def contactAcceptable(self, aminoacid, candidate):
        return candidate.symbol is Aminoacid.HYDROPHOBIC \
                and candidate.Id not in [aminoacid.Id-1, aminoacid.Id+1] 

    def calculateEnergy(self, eps=1.0):
        hContacts = 0
        for key in self.coords:
            aminoacid = self.coords.get(key)
            if aminoacid.symbol is Aminoacid.HYDROPHOBIC:
                candidates = [self.coords.get(tuple(np.add(key,x)),None) \
                        for x in self.directions]
                neighbors = [aa for aa in candidates \
                        if aa and self.contactAcceptable(aminoacid, aa)]
                hContacts += len(neighbors)
        self.setContacts(hContacts)
        self.energy = -eps*self.contacts

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
    
    def stat(self, energy, temp, kb):
        return math.exp(-energy/(kb*temp))

    def metropolis(self, chain, noSteps, state, temp, kb):
        microstates = []
        state.calculateEnergy()
        while noSteps:
            newState = state.transform()
            newState.calculateEnergy()
            aProbability = min(1, self.stat(newState.energy, \
                    temp, kb)/self.stat(state.energy, temp, kb))
            if random.random() < aProbability:
                state = newState
            microstates.append(state)
            noSteps -= 1
        return microstates


class Output(object):
    """ class responsible for plotting etc """

    @staticmethod
    def plotResult(x, y, name, labelx='', labely=''):
        plt.plot(x,y)
        plt.xlabel(labelx)
        plt.ylabel(labely)
        plt.savefig('%s.svg' % name)
        plt.close()


class Simulation(object):

    kb = 1.0
    
    def __init__(self, **kwargs):
        kwargs.update({'chain': Chain(kwargs.pop('sequence'))})
        self.__dict__.update(kwargs)
    
    def countCV(self, microstates, temp):
        energies = [x.energy for x in microstates]
        return np.var(energies)/(self.kb*(temp**2))

    def countAvI(self, microstates):
        """ get average moment of inertia """
        Is = []
        for state in microstates:
            cOfMass = np.mean(state.coords.keys(), axis=0)
            I = np.sum(np.power(np.subtract(state.coords.keys(),cOfMass),2),\
                    axis=0)
            Is.append(I)
        return np.sum(np.mean(Is), axis=0)             
    
    def plotResults(self, results, **kwargs): #TODO
        temp = kwargs.get('temp', [])
        if not temp:
            temp = [x/100.0 for x in xrange(int(self.tMax*100), \
                    int((self.tMin-self.delta)*100), int(-self.delta*100))]
        separated = zip(*results)
        Output.plotResult(temp, separated[0], 'cv' , 'T', 'Cv (T)')
        Output.plotResult(temp, separated[1], 'averageI' , 'T', '<I> (T)')
        for index, res in enumerate(separated[2]):
            a=dict((i,res.count(i)) for i in res)
            plt.bar(a.keys(),a.values())
            plt.xlabel('n of contacts')
            plt.ylabel('counts')
            plt.xlim(0, 10)
            plt.title('T=%s' % temp[index])
            plt.savefig('contacts_countT%s.png' % temp[index])
            plt.close()
    
    def saveTraj(self, result):
        traj = open("re_tMin.txt","wb")
        for step in result:
            for i, state in enumerate(step, 1):
                traj.write('MODEL%d\nCOMMENT:    nOfContacts=%d\n' % (i, \
                        state.contacts))
                lines=['']*len(self.chain)
                for key, value in state.coords.items():
                    lines[value.Id]=\
                            'ATOM      %d  CA  %s   %d  A     %r   %r   0.000\n'%\
                            (value.Id+1, \
                            'LYS' if value.symbol is Aminoacid.POLAR else 'ALA',
                            value.Id+1, float(key[0]), float(key[1]))
                traj.write(''.join(lines)+'ENDMDL\n')
        traj.close()

    def start(self):
        raise NotImplementedError


class ReplicaExchange(Simulation):
    
    #TODO tmax, tmin etc should not be saved as object attributes,
    # because they're no longer needed;
    # replicas could be objects
    def __init__(self, **kwargs):
        super(ReplicaExchange, self).__init__(**kwargs)
    
    def start(self):
        d = (self.tMax-self.tMin)/(self.replicas-1)
        temps = [self.tMin+d*x for x in xrange(self.replicas)]
        initialStates = \
                [Microstate(Microstate.getInitialCoords(self.chain))]*len(temps)
        length = self.length
        result = []
        while length:
            states = []
            for x, temp in enumerate(temps):
                states.append(Metropolis().metropolis(self.chain, 
                    self.noTransitions, initialStates[x], temp, self.kb))
            if not length % self.frequency:
                r = random.choice(range(self.replicas))
                i = [r-1, r] if r else [r, r+1] 
                replicas = states[i[0]:i[1]+1]
                delta = (1/(self.kb*temps[i[1]]) - 1/(self.kb*temps[i[0]]))\
                        *(states[i[0]][-1].energy - states[i[1]][-1].energy)
                prob = min(1, math.e**(-delta))
                if random.random() < prob:
                    states[i[0]] = replicas[1]
                    states[i[1]] = replicas[0]
                    print 'exchange!'
            initialStates = [state[-1] for state in states]
            result.append(states)
            length -= 1
        
        tMinResult=[x[0] for x in result]
        self.saveTraj(tMinResult)
        
        """ simulated annealing - like """
        lastResult=result[-1]
        cvs = [self.countCV(x, temps[y]) for y, x in enumerate(lastResult)]
        avis = [self.countAvI(x) for x in lastResult]
        contacts = [[n.contacts for n in m] for m in lastResult]
        self.plotResults(zip(cvs, avis, contacts), temp=temps)


class SimulatedAnnealing(Simulation):

    def __init__(self, **kwargs):
        super(SimulatedAnnealing, self).__init__(**kwargs)

    def start(self):
        result = []
        temp = self.tMax
        initialState = Microstate(Microstate.getInitialCoords(self.chain))
        while round(temp,2) >= self.tMin:
            print 'T = ', temp
            microstates = Metropolis().metropolis(self.chain, \
                    self.noTransitions, initialState, temp, self.kb)
            result.append((self.countCV(microstates, temp), \
                    self.countAvI(microstates), \
                    [n.contacts for n in microstates]))
            initialState = microstates[-1]
            temp -= self.delta
        print 'Plotting...'
        self.plotResults(result)


if __name__ == '__main__':
    main()
