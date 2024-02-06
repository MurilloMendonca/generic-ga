#include "GA.hpp"
#include <thread>
template <typename IndividualType,
          typename FitnessCalculationStrategy =
              IFitnessCalculationStrategy<IndividualType>,
          typename ViabilityStrategy = IViabilityStrategy<IndividualType>,
          typename IndividualGenerationStrategy =
              IIndividualGenerationStrategy<IndividualType>>
class ParallelGA : public GA<IndividualType, FitnessCalculationStrategy,
                             ViabilityStrategy, IndividualGenerationStrategy> {
private:
  int numThreads;

public:
  ParallelGA(int numThreads,
             FitnessCalculationStrategy fitnessCalculationStrategy,
             ViabilityStrategy viabilityStrategy,
             IndividualGenerationStrategy individualGenerationStrategy)
      : GA<IndividualType, FitnessCalculationStrategy, ViabilityStrategy,
           IndividualGenerationStrategy>(fitnessCalculationStrategy,
                                         viabilityStrategy,
                                         individualGenerationStrategy) {
    this->numThreads = numThreads;
  }
  ParallelGA(int numThreads,
             std::function<float(const IndividualType &)> fitnessCalculationStrategy,
     std::function<bool(const IndividualType &)> viabilityStrategy,
     std::function<IndividualType()> generationStrategy, int populationSize,
     float mutationRate, float elitismRate)
      : GA<IndividualType>(fitnessCalculationStrategy,
                                         viabilityStrategy,
                                         generationStrategy,
                                         populationSize, mutationRate,
                                         elitismRate) {
    this->numThreads = numThreads;
  }

  void parallelFitnessCalculation() {
    std::vector<std::thread> threads;
    int populationSize = this->population.size();
    int chunkSize = populationSize / numThreads;
    for (int i = 0; i < numThreads; i++) {
      int start = i * chunkSize;
      int end = (i + 1) * chunkSize;
      if (i == numThreads - 1) {
        end = populationSize;
      }
      threads.push_back(std::thread([this, start, end]() {
        for (int j = start; j < end; j++) {
        if (this->viabilityStrategy->isViable(this->population[j].getPhenotype())) {
          this->population[j].setFitness(
          this->fitnessCalculationStrategy->calculateFitness(this->population[j].getPhenotype())); 
        }
        else {
          this->population[j].setFitness(0);
        }
        }
      }));
    }
    for (int i = 0; i < numThreads; i++) {
      threads[i].join();
    }
  }
  void runWithThreads(int g) {
    for (int i = 0; i < g; i++) {
      std::vector<IndividualType> elite;
      for (int j = 0; j < this->populationSize * this->elitismRate; j++) {
        elite.push_back(this->population[j]);
      }
      this->crossoverPopulation();
      this->mutatePopulation();
      this->parallelFitnessCalculation();
      this->population.insert(this->population.end(), elite.begin(),
                              elite.end());
      std::sort(this->population.begin(), this->population.end());

      this->population.erase(this->population.begin() + this->populationSize,
                             this->population.end());
      // std::cout << "Generation " << i << " best fitness: " <<
      // population[0].getFitness() << std::endl;
    }
  }
};
