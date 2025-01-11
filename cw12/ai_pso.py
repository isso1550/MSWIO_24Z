import numpy as np
def pso(objective_function, num_particles=30, num_dimensions=2, num_iterations=100,
        inertia_weight=0.5, cognitive_coefficient=1.5, social_coefficient=1.5):
    
    # Inicjalizacja pozycji i prędkości cząstek
    particles_position = np.random.rand(num_particles, num_dimensions)
    particles_velocity = np.random.rand(num_particles, num_dimensions)

    # Inicjalizacja najlepszych pozycji cząstek i roju
    personal_best_position = particles_position.copy()
    personal_best_value = np.array([objective_function(p) for p in particles_position])
    global_best_position = personal_best_position[np.argmin(personal_best_value)]
    global_best_value = min(personal_best_value)

    # Inicjalizacja list do zbierania wyników
    sum_objective_values = []
    best_objective_values = []

    # Główna pętla algorytmu PSO
    for _ in range(num_iterations):
        iteration_sum = 0
        for i in range(num_particles):
            # Aktualizacja prędkości cząstki
            r1, r2 = np.random.rand(2)
            particles_velocity[i] = (inertia_weight * particles_velocity[i]
                                     + cognitive_coefficient * r1 * (personal_best_position[i] - particles_position[i])
                                     + social_coefficient * r2 * (global_best_position - particles_position[i]))

            # Aktualizacja pozycji cząstki
            particles_position[i] += particles_velocity[i]

            # Aktualizacja najlepszej pozycji cząstki i roju
            current_value = objective_function(particles_position[i])
            iteration_sum += current_value
            if current_value < personal_best_value[i]:
                personal_best_position[i] = particles_position[i]
                personal_best_value[i] = current_value
                if current_value < global_best_value:
                    global_best_position = particles_position[i]
                    global_best_value = current_value

        # Zbieranie wyników dla danej iteracji
        sum_objective_values.append(iteration_sum)
        best_objective_values.append(global_best_value)
    
    return sum_objective_values, best_objective_values