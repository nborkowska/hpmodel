#! /usr/bin/env python

import math
import numpy as np

class RotationMatrix(np.ndarray):

    def __new__(self, angle):
        return np.array([[round(math.cos(angle)), round(-math.sin(angle))],
            [round(math.sin(angle)), round(math.cos(angle))]])


class Rotation(object):

    def __init__(self, point):
        self.point = point
    
    def rotate(self):           #i to powinna byc metoda mikrostanu, no ale niekoniecznie bo to moze sobie obracac macierz po prostu, moze tu nie bedzie spr czy dozwolone 
        pass


class R90(Rotation):

    matrix = RotationMatrix(0.5*math.pi)


class R180(Rotation):

    matrix = RotationMatrix(math.pi)


class R270(Rotation):

    matrix = RotationMatrix(1.5*math.pi)


class Aminoacid(object):

    HYDROPHOBIC = 'H'
    POLAR = 'P'

    def __init__(self, symbol):
        self.symbol = symbol

    def isPolar(self):
        return self.symbol == POLAR


class Microstate(object):

    sequence = ''
    
    def __init__(self, coords, energy=10000):
        self.coords = coords
        self.energy = energy
    
    def setEnergy(self):
        pass

    def transform(self):
        pass
