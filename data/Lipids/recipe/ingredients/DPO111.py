#include as follow : execfile('pathto/DPO111.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
DPO111= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/DPO111.sph',
radii = [[3.9399999999999999, 1.52, 0.0, 3.2000000000000002, 1.8, 5.0800000000000001, 3.73, 4.3799999999999999]],
cutoff_boundary = 0,
Type = 'MultiSphere',
cutoff_surface = 0,
gradient = '',
jitterMax = [0.5, 0.5, 0.10000000000000001],
packingPriority = 0,
rotAxis = [0.0, 2.0, 1.0],
nbJitter = 5,
molarity = 1.0,
rotRange = 6.2831,
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/DPO111.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'DPO111',
positions = [[(-6.29, -2.9500000000000002, -10.130000000000001), (-0.93000000000000005, -0.78000000000000003, -20.960000000000001), (-1.1799999999999999, -1.6200000000000001, -23.030000000000001), (-9.5399999999999991, -2.04, -2.7200000000000002), (0.64000000000000001, 1.74, -22.940000000000001), (9.4600000000000009, 4.4800000000000004, -6.8300000000000001), (-1.3, -2.3300000000000001, -16.82), (4.2400000000000002, 0.92000000000000004, -15.42)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(DPO111)
