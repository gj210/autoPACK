#include as follow : execfile('pathto/LPO128.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
LPO128= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/recipe/ingredients/LPO128.sph',
radii = [[4.7300000000000004, 4.6699999999999999, 3.4700000000000002, 5.2199999999999998]],
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
meshFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/geoms/ingredients_1/LPO128.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'LPO128',
positions = [[(1.79, 1.3799999999999999, 19.449999999999999), (4.6399999999999997, 3.5099999999999998, 12.32), (-2.9300000000000002, -2.1499999999999999, 22.870000000000001), (3.1000000000000001, 8.6699999999999999, 6.7699999999999996)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(LPO128)
