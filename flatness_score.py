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