#include as follow : execfile('pathto/DPO112.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
DPO112= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/recipe/ingredients/DPO112.sph',
radii = [[5.5499999999999998, 4.96, 5.1600000000000001, 3.8599999999999999]],
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
meshFile = '/Users/ludo/DEV/autofill_googlesvn/data/Lipids/geoms/ingredients_1/DPO112.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'DPO112',
positions = [[(4.75, 2.3300000000000001, 3.6400000000000001), (-0.23999999999999999, -0.28999999999999998, 11.369999999999999), (-1.6899999999999999, 1.27, 19.149999999999999), (-0.67000000000000004, -6.0300000000000002, 3.9500000000000002)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(DPO112)
