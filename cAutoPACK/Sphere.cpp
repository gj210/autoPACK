/*
###############################################################################
#
# autoPACK Authors: Graham T. Johnson, Mostafa Al-Alusi, Ludovic Autin, Michel Sanner
#   Based on COFFEE Script developed by Graham Johnson between 2005 and 2010 
#   with assistance from Mostafa Al-Alusi in 2009 and periodic input 
#   from Arthur Olson's Molecular Graphics Lab
#
# autopack.cpp Authors: Ludovic Autin
#
# Translation from Python initiated March 15, 2010 by Ludovic Autin
#
#
# Copyright: Graham Johnson Ludovic Autin ©2010
#
# This file "autopack.cpp" is part of autoPACK, cellPACK.
#    
#    autoPACK is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    
#    autoPACK is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with autoPACK (See "CopyingGNUGPL" in the installation.
#    If not, see <http://www.gnu.org/licenses/>.
#
#
###############################################################################
@author: Graham Johnson, Ludovic Autin, & Michel Sanner
*/

#include "Sphere.h"

void sphere::setCount(float volume) {
    double nbr = this->molarity * volume * .000602;// #Mod by Graham 8/18/11
    int nbi = (int) nbr; //ceil, floor ?              #Mod by Graham 8/18/11
    double nbmod;    
    if (nbi == 0) 
        nbmod = nbr;
    else 
        nbmod = fmod(nbr, nbi);  //#Mod by Graham 8/18/11 ??
    double randval = rand();                 //#Mod by Graham 8/18/11
    if (nbmod >= randval)             //   #Mod by Graham 8/18/11
        nbi = (int)nbi+1;              //#Mod by Graham 8/18/11
    int nb = nbi;                        //#Mod by Graham 8/18/11

    std::cout << "#ingredient " << this->name << " nbMol " << this->nbMol << " nb "<< nb << " total " << nb + this->nbMol << std::endl;
    this->nbMol = nb + this->nbMol;
}

