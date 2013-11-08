# -*- coding: utf-8 -*-
"""
###############################################################################
#
# autoPACK Authors: Graham T. Johnson, Mostafa Al-Alusi, Ludovic Autin, Michel Sanner
#   Based on COFFEE Script developed by Graham Johnson between 2005 and 2010 
#   with assistance from Mostafa Al-Alusi in 2009 and periodic input 
#   from Arthur Olson's Molecular Graphics Lab
#
# ldSequence.py Authors: Ludovic Autin
###############################################################################
@author: Ludovic Autin, Graham Johnson,  & Michel Sanner
"""

import math

#==============================================================================
# First implementation came from QuasiRandom.py
# https://code.google.com/p/pistol/source/browse/trunk/Pistol/QuasiRandom.py?r=4
#==============================================================================
"""\
    QuasiRandom.py - Generators for quasirandom sequences
DESCRIPTION
    Use a Halton sequence to generate quasirandom sequences. 

Copyright (c) 2003 Richard P. Muller (rmuller@sandia.gov). All rights
reserved. See the LICENSE file for licensing details.
"""

def haltonterm(i,base=2):
    h=0
    fac =  1./base

    while i != 0:
        digit = i % base
        h = h + digit*fac
        i = (i-digit)/base
        fac = fac/base
    return h

def halton(terms,base=2):
    data = []
    for i in range(1,terms+1): data.append(haltonterm(i,base))
    return data

def halton2(terms,base1=2,base2=3):
    data = []
    for i in range(1,terms+1):
        data.append((haltonterm(i,base1),haltonterm(i,base2),0.0))
    return data

def halton3(terms,base1=2,base2=3,base3=5):
    data = []
    for i in range(1,terms+1):
        data.append((haltonterm(i,base1),haltonterm(i,base2),
                     haltonterm(i,base3)))
    return data

prime = [2.0,3.0,5.0,7.0,11.0,13.0,17.0,19.0,23.0,29.0,31.0,37.0,
         41.0,43.0,47.0,53.0,59.0,61.0,67.0,71.0,73.0,79.0,83.0,
         89.0,97.0,101.0,103.0,107.0,109.0,113.0,127.0,131.0,137.0,
         139.0,149.0,151.0,157.0,163.0,167.0,173.0]

class HaltonSequence:
    """Another algorithm to compute a d-dimensional Halton Sequence.
    This one gives the same results as the other."""
    def __init__(self,dim,atmost=10000):
        if dim > len(prime): raise "dim < %d " % len(prime)
        self.dim = dim

        self.err = 0.9/(atmost*prime[dim-1])

        self.prime = [0]*self.dim
        self.quasi = [0]*self.dim
        for i in range(self.dim):
            self.prime[i] = 1./prime[i]
            self.quasi[i] = self.prime[i]
        return

    def __call__(self):
        for i in range(self.dim):
            t = self.prime[i]
            f = 1.-self.quasi[i]
            g = 1.
            h = t
            while f-h < self.err:
                g = h
                h = h*t
            self.quasi[i] = g+h-f
        return self.quasi

class SobolSequence:
    "From the Gnu Scientific Library"
    def __init__(self,dim):
        self.max_dim = 40
        self.bit_count = 30
        self.primitive_polynomials = [
            1,     3,   7,  11,  13,  19,  25,  37,  59,  47,
            61,   55,  41,  67,  97,  91, 109, 103, 115, 131,
            193, 137, 145, 143, 241, 157, 185, 167, 229, 171,
            213, 191, 253, 203, 211, 239, 247, 285, 369, 299
            ]
        self.degree_table = [
            0, 1, 2, 3, 3, 4, 4, 5, 5, 5,
            5, 5, 5, 6, 6, 6, 6, 6, 6, 7,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 
            7, 7, 7, 7, 7, 7, 7, 8, 8, 8
            ]
        self.v_init = [
            [0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
             1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
             1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
             1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [0, 0, 1, 3, 1, 3, 1, 3, 3, 1,
             3, 1, 3, 1, 3, 1, 1, 3, 1, 3,
             1, 3, 1, 3, 3, 1, 3, 1, 3, 1,
             3, 1, 1, 3, 1, 3, 1, 3, 1, 3],
            [0, 0, 0, 7, 5, 1, 3, 3, 7, 5,
             5, 7, 7, 1, 3, 3, 7, 5, 1, 1,
             5, 3, 3, 1, 7, 5, 1, 3, 3, 7,
             5, 1, 1, 5, 7, 7, 5, 1, 3, 3],
            [0,  0,  0,  0,  0,  1,  7,  9, 13, 11,
             1,  3,  7,  9,  5, 13, 13, 11,  3, 15,
             5,  3, 15,  7,  9, 13,  9,  1, 11,  7,
             5, 15,  1, 15, 11,  5,  3,  1,  7,  9],
            [0,  0,  0,  0,  0,  0,  0,  9,  3, 27,
             15, 29, 21, 23, 19, 11, 25,  7, 13, 17,
             1, 25, 29,  3, 31, 11,  5, 23, 27, 19,
             21,  5,  1, 17, 13,  7, 15,  9, 31,  9],
            [0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
             0,  0,  0, 37, 33,  7,  5, 11, 39, 63,
             27, 17, 15, 23, 29,  3, 21, 13, 31, 25,
             9, 49, 33, 19, 29, 11, 19, 27, 15, 25],
            [0,   0,  0,  0,  0,  0,    0,  0,  0,   0,
             0,   0,  0,  0,  0,  0,    0,  0,  0,  13,
             33, 115, 41, 79, 17,  29, 119, 75, 73, 105,
             7,  59, 65, 21,  3, 113,  61, 89, 45, 107],
            [0, 0, 0, 0, 0, 0, 0, 0,  0,  0,
             0, 0, 0, 0, 0, 0, 0, 0,  0,  0,
             0, 0, 0, 0, 0, 0, 0, 0,  0,  0,
             0, 0, 0, 0, 0, 0, 0, 7, 23, 39]
            ] # end of self.v_init
        self.dim = dim
        #self.v_direction = [[0]*self.bit_count]*self.max_dim
        #self.v_direction = [[0]*self.max_dim]*self.bit_count #??
        self.v_direction = []
        for i in range(self.bit_count):
            array = []
            for j in range(self.max_dim):
                array.append(0)
            self.v_direction.append(array)
        

        # Initialize direction table in dimension 0.
        for k in range(self.bit_count): self.v_direction[k][0] = 1

        # Initialize in remaining dimensions.
        for i_dim in range(1,self.dim):
            poly_index = i_dim
            degree_i = self.degree_table[poly_index]
            includ = [0]*8
            
            # Expand the polynomial bit pattern to separate
            # components of the logical array includ[].
            p_i = self.primitive_polynomials[poly_index]
            for k in range(degree_i-1,-1,-1):
                includ[k] = ((p_i % 2) == 1)
                p_i /= 2

            # Leading elements for dimension i come from v_init[][].
            for j in range(degree_i):
                self.v_direction[j][i_dim] = self.v_init[j][i_dim]

            # Calculate remaining elements for this dimension,
            # as explained in Bratley+Fox, section 2.
            for j in range(degree_i, self.bit_count):
                newv = self.v_direction[j-degree_i][i_dim]
                ell = 1
                for k in range(degree_i):
                    ell *= 2
                    if includ[k]:
                        newv = newv ^ (ell*self.v_direction[j-k-1][i_dim])
                self.v_direction[j][i_dim] = newv
        # done with for i_dim in range(1,self.dim)

        # Multiply columns of v by appropriate power of 2
        ell = 1
        for j in range(self.bit_count-2,-1,-1):
            ell *= 2
            for i_dim in range(self.dim):
                self.v_direction[j][i_dim] = self.v_direction[j][i_dim]*ell
        # 1/(common denom of the elements in v_direction)

        self.last_denominator_inv = 1.0 /(2.0 * ell)

        # final setup
        self.sequence_count = 0
        self.last_numerator_vec = [0]*self.dim
        return

    def __call__(self):
        ell = 0
        c = self.sequence_count
        while 1:
            ell += 1
            if c % 2 == 1:
                c/=2
            else:
                break
        if ell > self.bit_count:
            raise "Sobol failed for %d" % self.sequence_count

        v = [0]*self.dim
        for i_dimension in range(self.dim):
            direction_i = self.v_direction[ell-1][i_dimension]
            old_numerator_i = self.last_numerator_vec[i_dimension]
            new_numerator_i = old_numerator_i ^ direction_i
            self.last_numerator_vec[i_dimension] = new_numerator_i
            v[i_dimension] = new_numerator_i*self.last_denominator_inv
        self.sequence_count += 1
        return v

kOneOverThree = (1.0 / 3.0)
kOneOverFive  = (1.0 / 5.0)


#==============================================================================
# Second implementation came from :
# http://www.cse.cuhk.edu.hk/~ttwong/papers/udpoint/udpoint.pdf
#==============================================================================
def SphereHammersley(n):
    p=0.0
    t=0.0
    st=0.0
    phi=0.0
    phirad=0.0
    k=0
    kk=0
    result=[]
#    pos=0
#    for (k=0, pos=0 ; k<n ; k++)
    for k in range(n):    
        t = 0
        p=0.5
        kk=k
#        for (p=0.5, kk=k ; kk ; p*=0.5, kk>>=1)
        while (kk):            
            if (kk & 1): #// kk mod 2 == 1
                t += p
            kk>>=1
            p*=0.5
        t = 2.0 * t - 1.0 #// map from [0,1] to [-1,1]
        phi = (k + 0.5) / n; #// a slight shift
        phirad = phi * 2.0 * math.pi; #// map to [0, 2 pi)
        st = math.sqrt(1.0-t*t);
        px = st * math.cos(phirad);
        py = st * math.sin(phirad);
        pz = t;
        result.append([px,py,pz])
    return result

def  SphereHalton(n, p2,marge=math.radians(45.0)):
    p=0.0
    t=0.0
    st=0.0
    phi=0.0
    phirad=0.0
    ip=0.0
    k=0
    kk=0
#    pos=0
    result=[]
    a=0
    for k in range(n):  
        #for (k=0, pos=0 ; k<n ; k++)
        t = 0;   
        p=0.5
        kk=k
#        for (p=0.5, kk=k ; kk ; p*=0.5, kk>>=1)
        while (kk):            
            if (kk & 1): #// kk mod 2 == 1
                t += p
            kk>>=1
            p*=0.5
        t = 2.0 * t - 1.0; #// map from [0,1] to [-1,1]
        st = math.sqrt(1.0-t*t);
        phi = 0;
        ip = 1.0/p2; #// inverse of p2
#        for (p=ip, kk=k ; kk ; p*=ip, kk/=p2) #// kk = (int)(kk/p2)
        p=ip
        kk=k        
        while(kk):
            a=kk % p2
            if (a):
                phi += a * p;
            kk/=p2
            p*=ip
        phirad = phi * 4.0 * marge; #// map from [0,0.5] to [0, 2 pi)
        px = st * math.cos(phirad);
        py = st * math.sin(phirad);
        pz = t;
        result.append([px,py,pz])
    return result

#see http://packinon.sourceforge.net/programs/index.html
#
def HaltonSequence( n,  b):
#/// return term i of the base b Halton sequence
#/// You can think of this as, for example, just a generalization of Heckbert's bit-
#/// reversal distribution trick.
#/// E.g., when b=3, write n as a base 3 number, digit 0 -> which third of interval the
#/// sample is in, 1 -> which third of that, 2 -> which third of that, etc.
#{
    result = 0
    ip = 1.0 / b    
    p = ip
    while (n > 0):
        result += (n % b) * p
        n = n / b
        p *= ip
    return result;

#==============================================================================
# third implementation come from fast object distribution article in Siggraph
# translate in python 
# http://www.andrewwillmott.com/s2007
#==============================================================================
class cHaltonSequence3:
#/// This calculates the Halton sequence incrementally
#/// for a base 2/3/5 triplet.
#//    cFXVector3 mPoint;  ///< Current sample point.
    def __init__(self):    
        self.mX=0.0
        self.mY=0.0
        self.mZ=0.0
        
        self.mBase2=0#uint32_t 
        self.mBase3=0
        self.mBase5=0
        
    def  inc(self,):
#    ///< Advance to next point in the sequence. Returns the index of this point. 
#    /////////////////////////////////////
#    // base 2       
        oldBase2 = self.mBase2 #uint32_t
        self.mBase2+=1;
        diff = self.mBase2 ^ oldBase2
    
#        // bottom bit always changes, higher bits
#        // change less frequently.
        s = 0.5;
    
#        // diff will be of the form 0*1+, i.e. one bits up until the last carry.
#        // expected iterations = 1 + 0.5 + 0.25 + ... = 2
        while True :
            if (oldBase2 & 1):
                self.mX -= s
            else:
                self.mX += s           
            s *= 0.5        
            diff = diff >> 1
            oldBase2 = oldBase2 >> 1
            if not diff :
                break    
        
#        /////////////////////////////////////
#        // base 3: use 2 bits for each base 3 digit.
        
        mask = 0x3  #// also the max base 3 digit
        add  = 0x1  #// amount to add to force carry once digit==3
        s = kOneOverThree
    
        self.mBase3+=1
    
#        // expected iterations: 1.5
        while True:
            if ((self.mBase3 & mask) == mask):
                self.mBase3 += add #         // force carry into next 2-bit digit
                self.mY -= 2 * s
                
                mask = mask << 2
                add  = add  << 2
                
                s *= kOneOverThree
            else :
                self.mY += s     #// we know digit n has gone from a to a + 1
                break
    
#        /////////////////////////////////////
#        // base 5: use 3 bits for each base 5 digit.
        mask = 0x7
        add  = 0x3#  // amount to add to force carry once digit==dmax
        dmax = 0x5#  // max digit
        
        s = kOneOverFive;
    
        self.mBase5+=1;
        
#        // expected iterations: 1.25
        while True:
            if ((self.mBase5 & mask) == dmax):
                self.mBase5 += add          #// force carry into next 3-bit digit
                self.mZ -= 4 * s
                
                mask = mask << 3
                dmax = dmax << 3
                add  = add  << 3
                
                s *= kOneOverFive
            else :
                self.mZ += s #     // we know digit n has gone from a to a + 1
                break
    
        return self.mBase2  #// return the index of this sequence point
    
    def reset(self,):
#    ///< Move back to first point in the sequence (i.e. the origin.)
        self.mBase2 = 0
        self.mBase3 = 0
        self.mBase5 = 0
        self.mX = 0.0
        self.mY = 0.0
        self.mZ = 0.0
