import numpy as np
import random

from .Chromosome import *

class Population:
    
    def __init__(self, popsize, n, given_matrix):
        
        self.elements = n
        self.popsize = popsize
        self.population = [Chromosome(self.elements, index=i) for i in range(popsize)]
        self.given_matrix = given_matrix
        
        
    def initializePopulation(self, function):

        """
        Initialize the population based on the fitness function selected to obtain the fitness value for each individual
        """

        for i, ind in enumerate(self.population):
            
            ind.randomInitialize(self.given_matrix, function, index = i)
         
        
    def printPopulation(self):

        """
        Print the __init__ variables for the population
        """

        print("Sudoku size: ", self.elements)
        print("Population size: ", self.popsize)

        for ind in self.population:
            print("\nIndividual: ")
            ind.printChromosome()
          
        
    def worst(self):
        
        """
        Obtain the worst individual, depending on the type_function
        """
        
        index_worst = max(range(len(self.population)), key=lambda i: self.population[i].fitness)
        
        return self.population[index_worst].fitness
    
    
    def best(self):
        
        """
        Obtain the best individual, depending on the type_function
        """

        index_best = min(range(len(self.population)), key=lambda i: self.population[i].fitness)

        return self.population[index_best].fitness
    
    
    # Obtain metrics the metrics of the population
    
    
    def totalFitness(self):

        """
        Obtain the total fitness value of the population
        """

        fitness_population = sum(ind.fitness for ind in self.population)

        return fitness_population

    def meanFitness(self):

        """
        Obtain the mean for the fitness of all members from the population
        """

        total_fit = self.totalFitness()

        mean_fitness = total_fit / self.popsize

        return mean_fitness

    def standarDeviation(self):

        """
        Obtain the standar deviation of the total fitness value in the population
        """

        fitness_population = [ind.fitness for ind in self.population]

        standarDev = np.std(fitness_population)

        return standarDev
    
    
#     def tournamentSelection(self, num_parents, tournament_size):
        
#         """
#         Select the parents to perform the crossover method
#         """
        
#         parents = []

#         for _ in range(num_parents):
            
#             selected_parent = self.selectUniqueParent(tournament_size, parents)
#             parents.append(selected_parent)

#         return parents
    

#     def selectUniqueParent(self, tournament_size, parents):
        
#         """
#         Select a parent and ensure that it is not in the parents list
#         """
        
#         flag = True
        
#         while flag:
            
#             # Shuffle the population for each tournament
            
#             random.shuffle(self.population)
            
#             # Perform the tournament
#             tournament = self.population[:tournament_size]
            

#             # Select the parent with the best fitness value
            
#             selected_parent = min(tournament, key=lambda ind: ind.fitness)
# #                 print(selected_parent.index)
            
#             if len(parents) > 0:
                
#                 # Check if the index is unique
                
# #                 print(selected_parent.index)
# #                 print(parents[0].index)
#                 if selected_parent.index != parents[0].index:
                
#                     flag = False

#                     return selected_parent
                
#             else:
                
#                 flag = False
                
#         return selected_parent
    
    
    def rows_crossover(self, parent1, parent2, cross_rate_rows):
        
        """
        Perform the crossover over the rows in the parents to obtain both childs
        """
        
        number_rows = self.elements
#         parent1, parent2 = parents[0], parents[1]
        
        # Create the two childs
        
        child1 = Chromosome(self.elements, 0, parent1.general_matrix.copy(), parent1.binary_matrix.copy())
        child2 = Chromosome(self.elements, 0, parent2.general_matrix.copy(), parent2.binary_matrix.copy())
        
        for i in range(number_rows):
            
            
            # Generate a random value to determine if the current row will be exchanged
            
            random_number = random.uniform(0.0, 1)
            
            if random_number <= cross_rate_rows:
                
                # Exchange the i-th row in both parents
                
                row_copy = child1.general_matrix[i].copy()
#                 print(row_copy)
                
                child1.general_matrix[i] = child2.general_matrix[i]
                child2.general_matrix[i] = row_copy
        
        return child1, child2
    
    
    def swap_mutation(self, mutation_rate_rows, mutation_rate_init):
        
        """
        Perform a mutation over one individual by swaping two positions inside random rows
        """
        
        
        for ind in self.population:
            
            # Iterate over each row from the individual
            
            for i in range(self.elements):
                
                if random.uniform(0.0, 1) < mutation_rate_rows:
                    
                    # Obtain the position of non-given numbers
                    
                    non_given_numbers = np.where(ind.binary_matrix[i] == 0)
#                     non_given_numbers = np.count_nonzero(ind.binary_matrix[i] == 0)
                    
                    if non_given_numbers[0].size >= 2:
                        
                        # Randomly select two different positions among non-given numbers
                        
                        index1, index2 = np.random.choice(non_given_numbers[0], size=2, replace=False)
                        
                        # Swap the values in the curent row
                        
                        aux_value = ind.general_matrix[i, index1]
                        
                        ind.general_matrix[i, index1] = ind.general_matrix[i, index2]
                        
                        ind.general_matrix[i, index2] = aux_value
                
                
                if random.uniform(0.0, 1) < mutation_rate_init:
                    
                    # Reinitialize the current row
                    
                    fixed_cols = np.where(self.given_matrix[i] != 0)[0]
            
                    # Generate values for the remaining values
                
                    remaining_values = np.setdiff1d(ind.values, self.given_matrix[i, fixed_cols])
                    shuffled_values = np.random.permutation(remaining_values)
            
                    # Assign the shuffled values to the corresponding columns
            
                    ind.general_matrix[i, self.given_matrix[i] == 0] = shuffled_values
                
                
    def update_population(self, mutation_rate_rows, mutation_rate_init, cross_rate, cross_rate_rows, function, tournament_size, elite_size):
 
        """
        Update population using tournament selection and mutation
        """     

        new_population = []
        elements = self.elements
        num_parents = 2

    
        # Crossover
        
        for ind in self.population:
            
            if random.uniform(0.0, 1) < cross_rate:
                
                # Select parents and perform the crossover method
                
                parent1 = ind
                parent2 = random.choice(self.population)
                
                child1, child2 = self.rows_crossover(parent1, parent2, cross_rate_rows)
                new_population.extend([child1, child2])

#         for _ in range(self.popsize // num_parents):
            
#             if random.uniform(0, 1) <= cross_rate:

#                 # Select parents and perform crossover using tournament selection
#                 parents = self.tournamentSelection(num_parents, tournament_size)

#                 parent1, parent2 = parents[0], parents[1]
                
#                 child1, child2 = self.rows_crossover(parent1, parent2, cross_rate_rows)
                                

#                 # Calculate the fitness value for each child
# #                 child1.evaluateFunction(function)
# #                 child2.evaluateFunction(function)

#                 new_population.extend([child1, child2])
            
        # Combine new children with the existing population
        combined_population = self.population + new_population
            
            
        # Update the current population
        self.population = combined_population
        
#         # Evaluate the new population
#         self.evaluate_population(function)
        
            
        # Mutation
        self.swap_mutation(mutation_rate_rows, mutation_rate_init)
        
        # Local search method
        
        self.local_search()
        
        # Evaluate the new population
        
        self.evaluate_population(function)
        

        # Sort the population based on the function_type
        sorted_population = sorted(self.population, key=lambda x: x.fitness, reverse=False)
        
        # Update the current population
        
        self.population = sorted_population[:self.popsize]


        # Update the index
        for i in range(0, len(self.population)):
            self.population[i].index = i
            
        # Elite population learning
        
        self.elite_learning(elite_size, function)
            
            
    def evaluate_population(self, function):
        """
        Obtain the fitness value for each individual in the population 
        """
        
        for ind in self.population:
            
            ind.evaluateFunction(function)
            
    # Column and Sub-Block Local Search
    
    def local_search(self):
        
        """
        Perform a local search over each individual in the population
        """
        
        for ind in self.population:
            
            ind.column_search()
        
        for ind in self.population:
            
            ind.sub_Block_search()
            
    # Elite learning
    
    def elite_learning(self, elite_size, function):
        
        """
        Define a elite population and replace the worst individual based on this population
        """
        
        # Define the elite population
        
        elite_population = self.population[:elite_size]
        
        # Obtain the probability Pb (probability of worst individual to be replaced for a individual on the elite_population)
        
        # Pb = (Maxfx - fx(xrandom) )/  Maxfx
        
        # Obtain the worst individual
        
        worst_individual = self.population[-1]
        
        # Select a randomm individual from the elite_population
        
        index_random = random.randint(0, elite_size-1)
        random_individual = elite_population[index_random]
        
#         random_individual.printChromosome()
#         worst_individual.printChromosome()
        
        # Calculate Pb
        
        Pb = (worst_individual.fitness - random_individual.fitness) / worst_individual.fitness
        
#         print("Pb", Pb)
        
        # Decide if the worst individual will be replaced or re-initialized
        
        if random.uniform(0.0, 1) < Pb:
            
            # Replace the worst individual
            
            print("Replace the worst individual", index_random)
            
            self.population[-1] = elite_population[index_random]
            
        else:
            
            # Initialize the worst individual
            
            print("Initialize the worst individual")
            self.population[-1].randomInitialize(self.given_matrix, function, index = worst_individual.index)