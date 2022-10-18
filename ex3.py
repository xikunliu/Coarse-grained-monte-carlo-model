#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 16:55:59 2019

@author: xigua
"""

import numpy as np
import time
import math
import random

#import compiled fortran library
import ex3lib
import sys
import os

# generate amino acid sequence AA  
# 1 = positive, 2 = negative, 3 = hydrophobic, 4 = hydrophilic,
# amino acid sequance for capped PHF6 is 3,3,4,3,3,3,1,3
capPHF6 = np.array([3,3,4,3,3,3,1,3])
capAA = np.tile (capPHF6,25)
# C terminal capped PHF6  
CAPHF6 = np.array([1,3,4,3,3,3,1,3])
CAAA = np.tile (CAPHF6,25)
# N terminal capped PHF6  
NAPHF6 = np.array([3,3,4,3,3,3,1,2])
AA = np.tile (CAPHF6,25)
# uncapped PHF6
PHF6 = np.array([1,3,4,3,3,3,1,2]) 
AAA = np.tile (PHF6,25)
#simulation for CA
#AA = capAA
#one atom size
sigma=2.5
#distance cutoff for pairwise interactions (angstrom) 2.5 sigma 
rc = 6.25
#set the random number seed; useful for debugging
np.random.seed(342324)
#number of amino acid per peptide
M = 8
#total numbe of amino acid
N = 200
#Length of the cubic box (angstrom)
L = 2.5*(N/0.8)**(1/3) 
#Temperature T
T = 300.0
#Bjerrum length BL (7 Å in water at room temperature)
BL = 7
#Screening length SL = 0.304/sqrt(Ionic strength 10 mM NaCl) 
SL = 30.4

def LineSearch(Pos, Dir, dx, EFracTol,
               Accel = 1.5, MaxInc = 10., MaxIter = 10000):
    """Performs a line search along direction Dir.
Input:
    Pos: starting positions, (N,3) array
    Dir: (N,3) array of gradient direction
    dx: initial step amount
    EFracTol: fractional energy tolerance
    Accel: acceleration factor
    MaxInc: the maximum increase in energy for bracketing
    MaxIter: maximum number of iteration steps
Output:
    PEnergy: value of potential energy at minimum along Dir
    Pos: minimum energy (N,3) position array along Dir
"""
    #identify global variables
    global M, L, rc, BL, SL, sigma 
    
    #start the iteration counter
    Iter = 0

    #find the normalized direction    
    NormDir = np.clip(Dir, -1.e100, 1.e100)
    NormDir = NormDir / np.sqrt(np.sum(NormDir * NormDir))

    #take the first two steps and compute energies    
    Dists = [0., dx]
    PEs = [ex3lib.calcenergy(Pos + NormDir * x, AA, M, L, rc, BL, SL, sigma) for x in Dists]
    
    #if the second point is not downhill in energy, back
    #off and take a shorter step until we find one
    while PEs[1] > PEs[0]:
        Iter += 1
        dx = dx * 0.5
        Dists[1] = dx
        PEs[1] = ex3lib.calcenergy(Pos + NormDir * dx, AA, M, L, rc, BL, SL, sigma)
        
    #find a third point
    Dists = Dists + [2. * dx]
    PEs = PEs + [ex3lib.calcenergy(Pos + NormDir * 2. * dx, AA, M, L, rc, BL, SL, sigma)]
    
    #keep stepping forward until the third point is higher
    #in energy; then we have bracketed a minimum
    while PEs[2] < PEs[1]:
        Iter += 1
            
        #find a fourth point and evaluate energy
        Dists = Dists + [Dists[-1] + dx]
        PEs = PEs + [ex3lib.calcenergy(Pos + NormDir * Dists[-1], AA, M, L, rc, BL, SL, sigma)]

        #check if we increased too much in energy; if so, back off
        if (PEs[3] - PEs[0]) > MaxInc * (PEs[0] - PEs[2]):
            PEs = PEs[:3]
            Dists = Dists[:3]
            dx = dx * 0.5
        else:
            #shift all of the points over
            PEs = PEs[-3:]
            Dists = Dists[-3:]
            dx = dx * Accel
            
    #we've bracketed a minimum; now we want to find it to high
    #accuracy
    OldPE3 = 1.e300
    while True:
        Iter += 1
        if Iter > MaxIter:
            print ("Warning: maximum number of iterations reached in line search.")
            break
            
        #store distances for ease of code-reading
        d0, d1, d2 = Dists
        PE0, PE1, PE2 = PEs

        #use a parobolic approximation to estimate the location
        #of the minimum
        d10 = d0 - d1
        d12 = d2 - d1
        Num = d12*d12*(PE0-PE1) - d10*d10*(PE2-PE1)
        Dem = d12*(PE0-PE1) - d10*(PE2-PE1)
        if Dem == 0:
            #parabolic extrapolation won't work; set new dist = 0 
            d3 = 0
        else:
            #location of parabolic minimum
            d3 = d1 + 0.5 * Num / Dem
            
        #compute the new potential energy
        PE3 = ex3lib.calcenergy(Pos + NormDir * d3, AA, M, L, rc, BL, SL, sigma)
        
        #sometimes the parabolic approximation can fail;
        #check if d3 is out of range < d0 or > d2 or the new energy is higher
        if d3 < d0 or d3 > d2 or PE3 > PE0 or PE3 > PE1 or PE3 > PE2:
            #instead, just compute the new distance by bisecting two
            #of the existing points along the line search
            if abs(d2 - d1) > abs(d0 - d1):
                d3 = 0.5 * (d2 + d1)
            else:
                d3 = 0.5 * (d0 + d1)
            PE3 = ex3lib.calcenergy(Pos + NormDir * d3, AA, M, L, rc, BL, SL, sigma)
            
        #decide which three points to keep; we want to keep
        #the three that are closest to the minimum
        if d3 < d1:
            if PE3 < PE1:
                #get rid of point 2
                Dists, PEs = [d0, d3, d1], [PE0, PE3, PE1]
            else:
                #get rid of point 0
                Dists, PEs = [d3, d1, d2], [PE3, PE1, PE2]
        else:
            if PE3 < PE1:
                #get rid of point 0
                Dists, PEs = [d1, d3, d2], [PE1, PE3, PE2]
            else:
                #get rid of point 2
                Dists, PEs = [d0, d1, d3], [PE0, PE1, PE3]
                
        #check how much we've changed
        if abs(OldPE3 - PE3) < EFracTol * abs(PE3):
            #the fractional change is less than the tolerance,
            #so we are done and can exit the loop
            break
        OldPE3 = PE3

    #return the position array at the minimum (point 1)        
    PosMin = Pos + NormDir * Dists[1]
    PEMin = PEs[1]
        
    return PEMin, PosMin

        
def ConjugateGradient(Pos, dx, EFracTolLS, EFracTolCG):
    """Performs a conjugate gradient search.
Input:
    Pos: starting positions, (N,3) array
    dx: initial step amount
    EFracTolLS: fractional energy tolerance for line search
    EFracTolCG: fractional energy tolerance for conjugate gradient
Output:
    PEnergy: value of potential energy at minimum
    Pos: minimum energy (N,3) position array 
"""
    #identify global variables
    global M, L, rc, BL, SL, sigma
    #initial search direction
    Forces = np.zeros_like(Pos)
    PE, Forces = ex3lib.calcenergyforces(Pos, AA, M, L, rc, BL, SL, sigma, Forces)
    Dir = Forces
    OldPE = 1.e300
    #iterative line searches
    while abs(PE - OldPE) > EFracTolCG * abs(PE):
        OldPE = PE
        PE, Pos = LineSearch(Pos, Dir, dx, EFracTolLS)
        OldForces = Forces.copy()
        PE, Forces = ex3lib.calcenergyforces(Pos, AA, M, L, rc, BL, SL, sigma, Forces)
        Gamma = np.sum((Forces - OldForces) * Forces) / np.sum(OldForces * OldForces)
        Dir = Forces + Gamma *  Dir
    return PE, Pos


def InitPositions(N, L):
    """Returns an array of initial positions of each atom,
placed on a cubic lattice for convenience.
Input:
    N: number of atoms
    L: box length
Output:
    Pos: (N,3) array of positions
"""
    #make the position array
    Pos = np.zeros((N,3), float)
    #compute integer grid # of locations for cubic lattice
    NLat = int(N**(1./3.) + 1.)
    LatSpac = L / NLat
    #make an array of lattice sites
    r = LatSpac * np.arange(NLat, dtype=float) - 0.5*L
    #loop through x, y, z positions in lattice until done
    #for every atom in the system
    i = 0
    for x in r:
        for y in r:
            for z in r:
                Pos[i] = np.array([x,y,z], float)
                #add a random offset to help initial minimization
                Offset = 0.1 * LatSpac * (np.random.rand(3) - 0.5)
                Pos[i] = Pos[i] + Offset
                i += 1
                #if done placing atoms, return
                if i >= N:
                    return Pos
    return Pos


def RescaleVelocities(Vel, T):
    """Rescales velocities in the system to the target temperature.
Input:
    Vel: (N,3) array of atomic velocities
    T: target temperature
Output:
    Vel: same as above
"""
    #recenter to zero net momentum (assuming all masses same)
    Vel = Vel - Vel.mean(axis=0)
    #find the total kinetic energy
    KE = 0.5 * np.sum(Vel * Vel)
    #find velocity scale factor from ratios of kinetic energy
    VScale = np.sqrt(1.5 * len(Vel) * T / KE)
    Vel = Vel * VScale
    return Vel  


def InitVelocities(N, T):
    """Returns an initial random velocity set.
Input:
    N: number of atoms
    T: target temperature
Output:
    Vel: (N,3) array of atomic velocities
"""
    Vel = np.random.rand(N, 3)
    Vel = RescaleVelocities(Vel, T)
    return Vel


def InitAccel(Pos):
    """Returns the initial acceleration array.
Input:
    Pos: (N,3) array of atomic positions
Output:
    Accel: (N,3) array of acceleration vectors
"""
    #global variables
    global M, L, rc, BL, SL, sigma
    #get the acceleration from the forces
    Forces = np.zeros_like(Pos)
    PEnergy, Accel = ex3lib.calcenergyforces(Pos, AA, M, L, rc, BL, SL, sigma, Forces)
    return Accel


def InstTemp(Vel):
    """Returns the instantaneous temperature.
Input:
    Vel: (N,3) array of atomic velocities
Output:
    Tinst: float
"""
    return np.sum(Vel * Vel) / (3. * len(Vel))

print('#simulation for CAPHF6')
#place atoms in a cubic coordinates 
Pos = InitPositions(N,L)    
# Minimization using the conjugate-gradiant method
NStepsMinimization=1000
for i in range(NStepsMinimization):
    PEnergy,Pos =ConjugateGradient(Pos, 0.001, 1.e-8, 1.e-10)
print(PEnergy)
# Equilibration using Monte Carlo simulation

def MCSimulation(Pos,PEnergy,MaxDisplacement):
    '''
    Returns energy and position after a monte carlo step
    Input:
        Pos:starting position, (N,3) array 
        PEnergy:starting potential energy
    Output:
        Pos: atomic position after a monte carlo step
        PEnergy:potential energy after a monte carlo step
        k:if the proposed position is accepted,k=1, else k=0
    '''
    # single atom displacement by a random amount
    MoveAtom = np.random.randint(N)
    OldPos = Pos[MoveAtom,:].copy()
    #compute old energy
    OldEnergy = ex3lib.calcenergyspecified(Pos, AA, M, L, rc, BL, SL, sigma, MoveAtom)
    #modify Pos[MoveAtom,:]
    Displacement=MaxDisplacement*L*np.random.uniform(low=-1, high=1, size=(3,))
    Pos[MoveAtom,:]=Pos[MoveAtom,:]+Displacement
    #compute new energy
    NewEnergy = ex3lib.calcenergyspecified(Pos, AA, M, L, rc, BL, SL, sigma, MoveAtom)
    #decide to accept or reject
    EnergyDiff=NewEnergy-OldEnergy
    if EnergyDiff > 0:
        Pacc=math.exp(-EnergyDiff/T)
        if Pacc <=random.randrange(0,1):
            Pos[MoveAtom,:]=OldPos
            k=0
        else:
            PEnergy=PEnergy+EnergyDiff
            k=1
    else:
        PEnergy = PEnergy+EnergyDiff
        k=1          
    return Pos,PEnergy,k

StepsEquil = 50000
fullUpdateStep = 10 * N
MaxDisplacement= 0.5 #max displacement that can achieve 80% exchange rate

print('#equilibration with change of number of contacts')
Cut = 5
#Contacts = np.zeros(N)


for i in range (StepsEquil):
    Pos,PEnergy,k = MCSimulation(Pos,PEnergy,MaxDisplacement)
    if i % fullUpdateStep == 0:
        PEnergy = ex3lib.calcenergy(Pos, AA, M, L, rc, BL, SL, sigma)
        ArrayContacts, Avg = ex3lib.calccontacts(Pos, L, Cut)
        with open ("ContactNumberENA.csv","a") as outfile:
            outfile.write(f"{i},{Avg}\n")
print('# Production')           
NStepsProd=50000
import atomwriteRevised
Names=["C" for i in range(N-M)]
for i in range(M):
    Names.append("O")
Pdb=atomwriteRevised.pdbfile("animNA.pdb",AtomNames=Names)

for i in range(NStepsProd):
    Pos,PEnergy,k = MCSimulation(Pos,PEnergy,MaxDisplacement)
    if i % fullUpdateStep == 0:
        PEnergy = ex3lib.calcenergy(Pos, AA, 25, L, rc, BL, SL, sigma)
        ArrayContacts, Avg = ex3lib.calccontacts(Pos, L, Cut)
        AvgsList=[]
        for NResidue in range (8):
            Avgs = ex3lib.calccontactsspecified(M, ArrayContacts, NResidue)
            AvgsList.append(Avgs)
        AvgsString=str(AvgsList).strip('[]')
        #calculate contact numbe of an atom of a specific type
        AvgAAList=[]
        for NAA in range (1,5):
       #     print(NAA)
       #     print("I am here")
            AvgAA = ex3lib.calccontactsaa(AA, ArrayContacts,NAA)
            AvgAAList.append(AvgAA)
       # print("I am here")
        AvgAAString=str(AvgAAList).strip('[]')
        with open ("ContactNumberPNA.csv","a") as outfile:
            outfile.write(f"{i},{Avg},{AvgsString},{AvgAAString}\n")
    if i % 200 == 0:  #update frequency =200
        # write pdb
        Pdb.write(Pos)
Pdb.close()