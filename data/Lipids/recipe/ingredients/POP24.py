#include as follow : execfile('pathto/POP24.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP24= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/recipe/ingredients/POP24.sph',
radii = [[7.3300000000000001, 7.1500000000000004, 5.0199999999999996, 3.79]],
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
meshFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/geoms/ingredients_1/POP24.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP24',
positions = [[(1.5900000000000001, 1.3500000000000001, 9.0700000000000003), (-1.1299999999999999, 1.3899999999999999, 20.649999999999999), (-8.1199999999999992, -2.4900000000000002, 7.46), (-3.9399999999999999, -1.1000000000000001, 16.140000000000001)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP24)
