#include as follow : execfile('pathto/LPO125.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
LPO125= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/recipe/ingredients/LPO125.sph',
radii = [[4.71, 8.2599999999999998, 5.79, 4.1299999999999999]],
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
meshFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/geoms/ingredients_1/LPO125.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'LPO125',
positions = [[(-1.78, 11.449999999999999, 11.300000000000001), (-4.1900000000000004, -5.8099999999999996, 8.9399999999999995), (-0.32000000000000001, -0.40999999999999998, 21.600000000000001), (-0.64000000000000001, 5.7699999999999996, 18.0)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(LPO125)
