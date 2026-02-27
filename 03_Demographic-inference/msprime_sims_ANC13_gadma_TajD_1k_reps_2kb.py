## Simulate summary stats based on best demographic model in GADMA
## Using msprime
## For ANC13 - North Central - best model
## I. Miranda, edited 2024-02-16


import msprime
import math

#Define best GADMA demographic model:

demography = msprime.Demography()

#Define model variables:
generation_time = 1
T_NC = 1007481 / generation_time  #time of split N-C

#Starting population sizes (current)
N_N = 367895      #present time pop size N
N_C = 2026672	  #present time pop size N

gr_N = -1.59677e-6   #grow rate from past to present - pop N
gr_C = 5.63597e-6   #grow rate from past to present - pop C


# Add all pops to the demographic model:

demography.add_population(
        name = "North",
        initial_size = N_N,
        growth_rate = gr_N,
        )

demography.add_population(
        name = "Central",
        initial_size = N_C,
        growth_rate = gr_C,
        )

demography.add_population (
        name = "Ancestral",
        initial_size = 693107,
        )

#Set initial migration events:

demography.set_migration_rate("North", "Central", 3.19e-7)  #set mig from N-C backwards in time, ie C-N forward in time
demography.set_migration_rate("Central", "North", 1.77e-6)  #set mig from C-N backwards in time, ie N-C forward in time


#Add Events that tell the model how populations correlate with each other:
demography.add_population_split(time = T_NC, derived = ["North", "Central"], ancestral = "Ancestral")


# Check demographic model definition:
print (demography.debug())



#Calculate summary stats from the data:

import numpy as np

#Define basic parameters:
n_reps = 1000
L_seq = 2000
mu = 2.2e-9


#Simulate n_rep trees for our demographic model:

def sim_replicates(num_replicates):
    ancestry_reps = msprime.sim_ancestry(
        samples = {"North": 111, "Central": 74},
        demography = demography,
        ploidy = 2,
        sequence_length = L_seq,
        random_seed = 12345,
        num_replicates = n_reps,
        )
    for ts in ancestry_reps:
        mutated_ts = msprime.sim_mutations(ts, rate = mu)
        yield mutated_ts



#Calculating Tajima's D North based on sites:

TajD_N_s = np.zeros(n_reps)

for replicate_index, ts in enumerate(sim_replicates(n_reps)):
        TajD_N_s[replicate_index] = ts.Tajimas_D (sample_sets = ts.samples(population=0))
TajD_N_s_mean = np.mean(TajD_N_s)
TajD_N_s_var = np.var (TajD_N_s)
print("Sites Tajima's D North mean:\t" + str(TajD_N_s_mean))
print("Sites Tajima's D North variance:\t" + str(TajD_N_s_var))

#Save Pi values to file:
file = open ("TajD_North_sites_sims_msprime_ANC13_best_gadma_1kreps.txt", "w+")

content = str(TajD_N_s)
file.write (content)
file.close()



#Calculating Tajima's D Central based on sites:

TajD_C_s = np.zeros(n_reps)

for replicate_index, ts in enumerate(sim_replicates(n_reps)):
        TajD_C_s[replicate_index] = ts.Tajimas_D (sample_sets = ts.samples(population=1))
TajD_C_s_mean = np.mean(TajD_C_s)
TajD_C_s_var = np.var (TajD_C_s)
print("Sites Tajima's D Central mean:\t" + str(TajD_C_s_mean))
print("Sites Tajima's D Central variance:\t" + str(TajD_C_s_var))

#Save Pi values to file:
file = open ("TajD_Central_sites_sims_msprime_ANC13_best_gadma_1kreps.txt", "w+")

content = str(TajD_C_s)
file.write (content)
file.close()
