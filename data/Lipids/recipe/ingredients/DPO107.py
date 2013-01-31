#include as follow : execfile('pathto/DPO107.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
DPO107= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/recipe/ingredients/DPO107.sph',
radii = [[4.7800000000000002, 2.3300000000000001, 4.1500000000000004, 8.0399999999999991]],
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
meshFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/geoms/ingredients_1/DPO107.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'DPO107',
positions = [[(-2.5800000000000001, 0.54000000000000004, 16.239999999999998), (-1.28, 2.1099999999999999, 27.109999999999999), (0.11, 0.90000000000000002, 22.510000000000002), (-3.4500000000000002, 0.94999999999999996, 7.29)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(DPO107)
