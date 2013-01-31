#include as follow : execfile('pathto/POP87.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP87= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/recipe/ingredients/POP87.sph',
radii = [[5.9199999999999999, 5.7000000000000002, 6.1699999999999999, 5.5300000000000002]],
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
meshFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/geoms/ingredients_1/POP87.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP87',
positions = [[(-0.17999999999999999, 3.6299999999999999, 21.300000000000001), (-1.5700000000000001, -0.72999999999999998, 23.390000000000001), (-5.5, -5.9800000000000004, 12.92), (-0.76000000000000001, 5.0599999999999996, 8.9800000000000004)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP87)
