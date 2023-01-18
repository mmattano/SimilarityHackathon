
'''

This file provides function for random pertrubation of spectral peaks

'''


import random
import numpy as np
from matchms import Spectrum


def normal_weights(x, N):
    return np.exp(-0.5 * (x / (N/3))**2) + 0.5

def draw_perturbation_signs(N=6):
    '''
    Producing an N-permutation of set {-1, 1}. 
    The difference between positiv and negative elements in terms of absolut numbers
    follow a normal distribution to only favor similar numbers of opposite signs (small differences). 
    
    Note that the difference between positiv and negative elements % 2 = N % 2 (as when one sign flips the difference increase by 2)
    
    '''
    candidates = [x for x in range(N + 1) if x % 2 == N%2] 
    weights = [normal_weights(d, N) for d in candidates]
    
    # Draw desired difference between the number of positive and negative signs
    draw = random.choices(candidates, k=1,
                weights=weights)

    # Select signs s.t. difference is fulfilled
    num_ones = int((N + draw[0]) / 2)
    signs = [1] * num_ones + [-1] * (N - num_ones)

    # Randomize order of signs 
    random.shuffle(signs)

    # multiply by +/- 1
    return random.choices([-1, 1], k=1)[0] * np.array(signs)



def apply_rand_perturbation(N, values, signs, perturb_range=[0, 0.5]):
    '''
    Perturb values by (1 +- r), where r is the perturbation strength randomly sampled from provided range (perturb_range).
    Signs need to be provided to determine the direction of the perturbation.
    '''
    
    # Draw perturbation strength r from uniform distbution (for all values)
    r = np.random.uniform(low=perturb_range[0], high=perturb_range[1], size=N)
    
    # Factor is (1 +/- r)
    perturbation_factors = 1 + signs * r 

    # Perturb values
    return values * perturbation_factors


def perturb_spetrum_copy(spectrum, perturb_range=[0, 0.5]):
    '''
    Creates a copy of the input spectrum with perturbed intensity values
    according to the perturbation range.
    '''
    
    new_intensities = apply_rand_perturbation(N=len(spectrum.intensities), values=spectrum.intensities, signs=draw_perturbation_signs(N=len(spectrum.intensities)), perturb_range=perturb_range)

    new_spectrum = Spectrum(mz=spectrum.mz,
                    intensities=new_intensities,
                    metadata={'id': 'None', 'comment': 'this is a copy', 'precursor_mz': spectrum.metadata["precursor_mz"]})

    return new_spectrum




