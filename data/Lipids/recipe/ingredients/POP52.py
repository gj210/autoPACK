#include as follow : execfile('pathto/POP52.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP52= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/recipe/ingredients/POP52.sph',
radii = [[5.9400000000000004, 6.0599999999999996, 4.6299999999999999, 4.4400000000000004]],
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
meshFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/geoms/ingredients_1/POP52.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'POP52',
positions = [[(7.3200000000000003, -4.1299999999999999, -6.9299999999999997), (1.6499999999999999, -0.79000000000000004, -18.289999999999999), (-4.8300000000000001, -1.72, -14.949999999999999), (-4.4500000000000002, -5.5899999999999999, -6.0599999999999996)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP52)
