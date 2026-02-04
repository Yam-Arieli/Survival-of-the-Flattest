import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def calculate_flatness(fitness, nei_fitness):
    """
    Calculates the flatness (robustness) score given a sequence's fitness 
    and the fitness values of its neighbors.
    
    Args:
        fitness (float): The fitness value of the wild-type sequence.
        nei_fitness (list[float]): A list of fitness values for all neighbors.
        
    Returns:
        float: The flatness score. 
               Lower is flatter (more robust).
               Higher is pointier (more fragile).
    """
    deleterious_drops = []
    
    for val in nei_fitness:
        # Calculate the drop in fitness
        drop = fitness - val
        
        # Apply "Rectified" logic:
        # If drop > 0 (neighbor is worse), we count it as a risk.
        # If drop <= 0 (neighbor is better), we ignore it (set to 0).
        if drop > 0:
            deleterious_drops.append(drop)
        else:
            deleterious_drops.append(0.0)
            
    # The score is the average risk of a random mutation
    flatness_score = sum(deleterious_drops) / len(nei_fitness)
    
    return flatness_score

def get_hamming_distance(s1, s2):
    """Returns the number of positions where s1 and s2 differ."""
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def get_neighbors(seq_list):
    """
    Takes a list of DNA sequences.
    For each sequence, returns the original followed by all 240 single-mutant neighbors.
    """
    bases = ['A', 'C', 'G', 'T']
    all_neighbors = []
    
    for group_idx, parent_seq in enumerate(seq_list):
        # Optional: Print a header to know which group this is
        all_neighbors.append(parent_seq)
        
        # Generate and print all 240 neighbors
        neighbor_count = 0
        for i, char in enumerate(parent_seq):
            for base in bases:
                if base != char:
                    # Create mutant (Hamming dist 1)
                    mutant = parent_seq[:i] + base + parent_seq[i+1:]
                    all_neighbors.append(mutant)
                    neighbor_count += 1
    
    return all_neighbors

def generate_sequences(n, length=80, min_distance=30):
    """
    Generates n distinct DNA sequences of a given length.
    Ensures that every new sequence is at least 'min_distance' 
    away from all previously generated ones.
    """
    bases = ['A', 'C', 'G', 'T']
    sequences = []
    attempts = 0
    max_attempts = n * 100  # Avoid infinite loops
    
    while len(sequences) < n and attempts < max_attempts:
        # 1. Generate a random candidate
        candidate = "".join(random.choices(bases, k=length))
        
        # 2. Check distance against all existing sequences
        is_far_enough = True
        for existing in sequences:
            if get_hamming_distance(candidate, existing) < min_distance:
                is_far_enough = False
                break
        
        # 3. Add or discard
        if is_far_enough:
            sequences.append(candidate)
        
        attempts += 1

    if len(sequences) < n:
        print(f"Warning: Could only find {len(sequences)} sequences with that distance constraint.")
        
    return sequences

def get_fit_from_interpolation(df_exp: pd.DataFrame, df_curve: pd.DataFrame):
    # 1. Sort the curve data (Crucial: interpolation requires sorted X values)
    df_curve = df_curve.sort_values('expression')

    # 2. Interpolate
    # np.interp(x_to_predict, x_reference, y_reference)
    return np.interp(
        df_exp['expression'],              # The values you have
        df_curve['expression'],        # The X axis of the curve
        df_curve['fitness']            # The Y axis of the curve
    )

def get_fit_n_flattness(df_exp: pd.DataFrame, seq_num: int):
    temp_df = df_exp[df_exp['group_id'] == seq_num]
    this_seq_fit = temp_df['fitness'].iloc[0]
    neighbors_fits = temp_df['fitness'].iloc[1:].tolist()
    return this_seq_fit, calculate_flatness(this_seq_fit, neighbors_fits)

def plot_fitness_vs_flatness(fit_flatt_list, power=4):
    fitness_values = [item[0]**power for item in fit_flatt_list]
    flatness_values = [item[1] for item in fit_flatt_list]
    first_seq = fit_flatt_list[0]

    plt.scatter(fitness_values, flatness_values, c='blue', alpha=0.5, label='all sequences')
    plt.scatter([first_seq[0]**power], [first_seq[1]], c='red', alpha=1, label='first sequence')
    plt.xlabel('Fitness')
    plt.ylabel('Flatness')
    plt.title('Fitness vs Flatness')
    plt.legend()
    plt.show()