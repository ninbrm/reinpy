## Project plan

### Habitat functionality

**GOAL:** 
Write paper for an ecological journal "Integrating the effects of habitat loss and fragmentation through graph-based habitat functionality" (working title). The idea is to illustrate the behavior of the model on a set of synthetic demonstration landscapes. We want to highlight the effects of habitat loss, fragmentation and their interaction; and how this is captured in habitat functionality. In the second part of the paper we will use the real case of Snøhetta (probably with the data and model from the JAE paper) to demonstrate the improved computation of the affinities, the computational performance (and tricks to improve it, e.g. reduce edges), and its applicability to real-world cases. 

**STEPS:**

1. Finish affinity script (using GRASS)
2. Create artificial landscapes: using the sum of two bivariate normal distributions and by varying the volume, mean, variance and covariance of each distribution we can generate a whole range of fragmentation scenarios. We can go from a circle to a oblong and star shaped range. In addition, we can go from two connected range to two disconnected ranges (and vary the relative size of each range). Alternatively we can use the ```RandomFields``` package to generate landscapes with different characteristics (see Chloe et al.). There we will have two/three parameters that we can control.
To simulate the model, the continuous variation in habitat quality is discretized onto a grid. Landscapes were generated with a Gaussian random field (GRF) using the RandomFields R package [38] to be 50 × 50 in size using an exponential covariance function with variance = 1, nugget = 0, and a set of mean (patch concentration, μQ) and scale (patch size, γQ) values (Fig 1, Table 1). Differing levels of patchiness were simulated by truncating all negative values to zero, resulting in more areas of no resources for smaller mean values. Landscapes were normalized to sum to one so that the total resources are the same across all landscapes. For each combination of μQ and γQ, twenty landscape replicates were generated. 


3. blabla
