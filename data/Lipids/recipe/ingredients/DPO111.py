#include as follow : execfile('pathto/DPO111.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
DPO111= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/recipe/ingredients/DPO111.sph',
radii = [[6.3600000000000003, 5.9100000000000001, 7.2999999999999998, 2.9399999999999999]],
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
meshFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/geoms/ingredients_1/DPO111.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'DPO111',
positions = [[(9.1500000000000004, 6.4900000000000002, -7.8200000000000003), (0.72999999999999998, 1.1299999999999999, -16.640000000000001), (-7.6299999999999999, -0.059999999999999998, -6.4900000000000002), (0.46000000000000002, 3.3300000000000001, -22.559999999999999)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(DPO111)
