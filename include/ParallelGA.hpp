#include "GA.hpp"
#include <thread>
#include <boost/asio.hpp>
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
  std::unique_ptr<boost::asio::thread_pool> pool;

public:
  ParallelGA(int numThreads,
             std::unique_ptr<FitnessCalculationStrategy> fitnessCalculationStrategy,
             std::unique_ptr<ViabilityStrategy> viabilityStrategy,
             std::unique_ptr<IndividualGenerationStrategy> individualGenerationStrategy,
             int populationSize, float mutationRate, float elitismRate)
      : GA<IndividualType, FitnessCalculationStrategy, ViabilityStrategy,
           IndividualGenerationStrategy>(std::move(fitnessCalculationStrategy),
                                       std::move(viabilityStrategy),
                                       std::move(individualGenerationStrategy),
                                       populationSize, mutationRate,
                                       elitismRate) {
    this->numThreads = numThreads;
    pool = std::make_unique<boost::asio::thread_pool>(numThreads);
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
    pool = std::make_unique<boost::asio::thread_pool>(numThreads);
  }

  void parallelFitnessCalculation() {
    int populationSize = this->population.size();
    int chunkSize = populationSize / numThreads;
    for (int i = 0; i < numThreads; i++) {
      int start = i * chunkSize;
      int end = (i + 1) * chunkSize;
      if (i == numThreads - 1) {
        end = populationSize;
      }
      boost::asio::post(*pool,[this, start, end]() {
        for (int j = start; j < end; j++) {
        if (this->viabilityStrategy->isViable(this->population[j].getPhenotype())) {
          this->population[j].setFitness(
          this->fitnessCalculationStrategy->calculateFitness(this->population[j].getPhenotype())); 
        }
        else {
          this->population[j].setFitness(0);
        }
        }
      });
    }
    pool->join();
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
