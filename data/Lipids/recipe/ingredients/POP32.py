#include as follow : execfile('pathto/POP32.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP32= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/recipe/ingredients/POP32.sph',
radii = [[2.8599999999999999, 5.7599999999999998, 5.8499999999999996, 5.1900000000000004]],
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
meshFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/geoms/ingredients_1/POP32.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP32',
positions = [[(2.6499999999999999, 2.0699999999999998, -21.379999999999999), (-0.34999999999999998, -0.79000000000000004, -16.859999999999999), (-0.96999999999999997, 5.3399999999999999, -6.2599999999999998), (3.9199999999999999, -3.71, -6.6399999999999997)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP32)
