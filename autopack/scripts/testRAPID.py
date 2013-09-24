# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 20:25:27 2013

@author: ludo
"""

#test rapid colllisiob detecion in autopack
import c4d
h=c4d.af.values()[0].histoVol.values()[0]
ingr=h.organelles[0].surfaceRecipe.ingredients[0]
h.rIngr