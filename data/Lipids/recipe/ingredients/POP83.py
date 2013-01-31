#include as follow : execfile('pathto/POP83.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP83= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/recipe/ingredients/POP83.sph',
radii = [[4.6500000000000004, 4.1399999999999997, 4.2999999999999998, 4.0300000000000002]],
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
meshFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/geoms/ingredients_1/POP83.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP83',
positions = [[(-0.70999999999999996, 2.6000000000000001, -16.140000000000001), (-0.20000000000000001, -2.8199999999999998, -9.7699999999999996), (3.1099999999999999, -1.96, -17.629999999999999), (4.6600000000000001, -5.7599999999999998, -5.5899999999999999)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP83)
