import numpy as np

def harmony_search(objective_function, num_harmonies=30, num_dimensions=2, num_iterations=100,
                   harmony_memory_consideration_rate=0.9, pitch_adjustment_rate=0.3, bandwidth=0.01):
    """
    Performs Harmony Search (HS) on the given objective function.

    Args:
        objective_function (callable): The objective function to be minimized.
        num_harmonies (int): The number of harmonies in the harmony memory. Default is 30.
        num_dimensions (int): The number of dimensions in the search space. Default is 2.
        num_iterations (int): The number of iterations to run the HS algorithm. Default is 100.
        harmony_memory_consideration_rate (float): The rate of considering the harmony memory. Default is 0.9.
        pitch_adjustment_rate (float): The rate of pitch adjustment. Default is 0.3.
        bandwidth (float): The bandwidth for pitch adjustment. Default is 0.01.

    Returns:
        list: Sum of objective values for each iteration.
        list: Best objective values for each iteration.
    """
    # Inicjalizacja pamięci harmonii
    harmony_memory = np.random.rand(num_harmonies, num_dimensions)
    harmony_memory_values = np.array([objective_function(harmony) for harmony in harmony_memory])

    # Inicjalizacja list do zbierania wyników
    sum_objective_values = []
    best_objective_values = []

    # Główna pętla algorytmu HS
    for _ in range(num_iterations):
        new_harmony = np.zeros(num_dimensions)
        for i in range(num_dimensions):
            if np.random.rand() < harmony_memory_consideration_rate:
                new_harmony[i] = harmony_memory[np.random.randint(num_harmonies), i]
                if np.random.rand() < pitch_adjustment_rate:
                    new_harmony[i] += bandwidth * (np.random.rand() - 0.5)
            else:
                new_harmony[i] = np.random.rand()
        
        new_harmony_value = objective_function(new_harmony)
        worst_index = np.argmax(harmony_memory_values)
        if new_harmony_value < harmony_memory_values[worst_index]:
            harmony_memory[worst_index] = new_harmony
            harmony_memory_values[worst_index] = new_harmony_value
        
        # Zbieranie wyników dla danej iteracji
        sum_objective_values.append(np.sum(harmony_memory_values))
        best_objective_values.append(np.min(harmony_memory_values))
    
    return sum_objective_values, best_objective_values