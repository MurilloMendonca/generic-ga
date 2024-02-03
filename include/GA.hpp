#pragma once
#include <algorithm>
#include <functional>
#include <memory>
#include <random>
#include <vector>
template <typename IndividualType> class IFitnessCalculationStrategy {
public:
  virtual ~IFitnessCalculationStrategy() = default;
  virtual float calculateFitness(const IndividualType &individual) = 0;
  virtual float operator()(const IndividualType &individual) {
    return calculateFitness(individual);
  }
};

template <typename IndividualType> class IViabilityStrategy {
public:
  virtual ~IViabilityStrategy() = default;
  virtual bool isViable(const IndividualType &individual) = 0;
  virtual bool operator()(const IndividualType &individual) {
    return isViable(individual);
  }
};

template <typename IndividualType> class IIndividualGenerationStrategy {
public:
  virtual ~IIndividualGenerationStrategy() = default;
  virtual IndividualType generateIndividual() = 0;
  virtual IndividualType operator()() { return generateIndividual(); }
};

template <typename IndividualType,
          typename FitnessCalculationStrategy =
              IFitnessCalculationStrategy<IndividualType>,
          typename ViabilityStrategy = IViabilityStrategy<IndividualType>,
          typename IndividualGenerationStrategy =
              IIndividualGenerationStrategy<IndividualType>>
class GA {
private:
  std::vector<IndividualType> population;
  float mutationRate;
  float elitismRate;
  int populationSize;

  std::unique_ptr<FitnessCalculationStrategy> fitnessCalculationStrategy;
  std::unique_ptr<ViabilityStrategy> viabilityStrategy;
  std::unique_ptr<IndividualGenerationStrategy> individualGenerationStrategy;
  void initializePopulation() {
    for (int i = 0; i < populationSize; i++) {
      population.push_back(
          this->individualGenerationStrategy->generateIndividual());
    }
    evaluatePopulation();
    std::sort(population.begin(), population.end());
  }
  void evaluatePopulation() {
    for (auto &individual : population) {
      if (!viabilityStrategy->isViable(individual.getPhenotype())) {
        individual.setFitness(0);
        continue;
      }
      individual.setFitness(fitnessCalculationStrategy->calculateFitness(
          individual.getPhenotype()));
    }
  }

  void mutatePopulation() {
    for (auto &individual : population) {
      mutation(individual);
    }
  }
  void crossoverPopulation() {
    std::vector<IndividualType> newPopulation;
    while (newPopulation.size() < populationSize) {
      std::pair<IndividualType, IndividualType> parents = selection();
      std::pair<IndividualType, IndividualType> children =
          crossover(parents.first, parents.second);
      newPopulation.push_back(children.first);
      newPopulation.push_back(children.second);
    }
    population = newPopulation;
  }

  std::pair<IndividualType, IndividualType> selection() {
    // Get a pair of distinct Individuals using tournament selection
    std::pair<IndividualType, IndividualType> parents;
    std::vector<IndividualType> tournament;
    // set up random seed
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, populationSize - 1);

    tournament.clear();
    for (int j = 0; j < 10; j++) {
      tournament.push_back(population[dis(gen)]);
    }
    std::sort(tournament.begin(), tournament.end());
    parents.first = tournament[0];
    parents.second = tournament[1];
    int i = 2;
    while (parents.first == parents.second && i < 10) {
      parents.second = tournament[i++];
    }

    return parents;
  }
  std::pair<IndividualType, IndividualType> crossover(IndividualType ind1,
                                                      IndividualType ind2) {
    std::vector<bool> genes1 = ind1.getGenotype();
    std::vector<bool> genes2 = ind2.getGenotype();
    // set up random seed
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_int_distribution<> dis(0, genes1.size() - 1);
    int crossoverPoint = dis(gen);

    std::vector<bool> newGenes1(genes1.begin(), genes1.begin() + crossoverPoint);
    newGenes1.insert(newGenes1.end(), genes2.begin() + crossoverPoint,
                     genes2.end()); 

    std::vector<bool> newGenes2(genes2.begin(), genes2.begin() + crossoverPoint);
    newGenes2.insert(newGenes2.end(), genes1.begin() + crossoverPoint,
                     genes1.end());

    return std::make_pair(IndividualType(newGenes1), IndividualType(newGenes2));
  }
  void mutation(IndividualType &i) {
    std::vector<bool> genes = i.getGenotype();
    // set up random seed
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> mut(0, 1);

    for (int j = 0; j < genes.size(); j++) {
      if (mut(gen) < mutationRate) {
        genes[j] = !genes[j];
      }
    }

    i.setGenotypeAndUpdatePhenotype(genes);
  }

public:
  // Compile-time strategy constructor
  GA(int populationSize, float mutationRate, float elitismRate)
      : populationSize(populationSize), mutationRate(mutationRate),
        elitismRate(elitismRate) {

    this->fitnessCalculationStrategy =
        std::make_unique<FitnessCalculationStrategy>();
    this->viabilityStrategy = std::make_unique<ViabilityStrategy>();
    this->individualGenerationStrategy =
        std::make_unique<IndividualGenerationStrategy>();
    initializePopulation();
  }
  GA() {
    this->populationSize = 100;
    this->mutationRate = 0.01;
    this->elitismRate = 0.1;
    this->fitnessCalculationStrategy =
        std::make_unique<FitnessCalculationStrategy>();
    this->viabilityStrategy = std::make_unique<ViabilityStrategy>();
    this->individualGenerationStrategy =
        std::make_unique<IndividualGenerationStrategy>();
    initializePopulation();
  }

  // Runtime strategy constructor
  GA(std::unique_ptr<FitnessCalculationStrategy> fitnessStrategy,
     std::unique_ptr<ViabilityStrategy> viabilityStrategy,
     std::unique_ptr<IndividualGenerationStrategy> generationStrategy,
     int populationSize, float mutationRate, float elitismRate)
      : fitnessCalculationStrategy(std::move(fitnessStrategy)),
        viabilityStrategy(std::move(viabilityStrategy)),
        individualGenerationStrategy(std::move(generationStrategy)),
        populationSize(populationSize), mutationRate(mutationRate),
        elitismRate(elitismRate) {
    initializePopulation();
  }

  GA(std::unique_ptr<FitnessCalculationStrategy> fitnessStrategy,
     std::unique_ptr<ViabilityStrategy> viabilityStrategy,
     std::unique_ptr<IndividualGenerationStrategy> generationStrategy)
      : fitnessCalculationStrategy(std::move(fitnessStrategy)),
        viabilityStrategy(std::move(viabilityStrategy)),
        individualGenerationStrategy(std::move(generationStrategy)) {
    this->populationSize = 100;
    this->mutationRate = 0.01;
    this->elitismRate = 0.1;
    initializePopulation();
  }

  // function strategy constructor
  GA(std::function<float(const IndividualType &)> fitnessCalculationStrategy,
     std::function<bool(const IndividualType &)> viabilityStrategy,
     std::function<IndividualType()> generationStrategy, int populationSize,
     float mutationRate, float elitismRate)
      : populationSize(populationSize), mutationRate(mutationRate),
        elitismRate(elitismRate) {

    class ConcreteFitnessCalculationStrategy
        : public IFitnessCalculationStrategy<IndividualType> {
      std::function<float(const IndividualType &)> fitnessCalculationStrategy;

    public:
      ConcreteFitnessCalculationStrategy(
          std::function<float(const IndividualType &)>
              fitnessCalculationStrategy)
          : fitnessCalculationStrategy(fitnessCalculationStrategy) {}
      float calculateFitness(const IndividualType &individual) override {
        return fitnessCalculationStrategy(individual);
      }
    };

    class ConcreteViabilityStrategy
        : public IViabilityStrategy<IndividualType> {
      std::function<bool(const IndividualType &)> viabilityStrategy;

    public:
      ConcreteViabilityStrategy(
          std::function<bool(const IndividualType &)> viabilityStrategy)
          : viabilityStrategy(viabilityStrategy) {}
      bool isViable(const IndividualType &individual) override {
        return viabilityStrategy(individual);
      }
    };

    class ConcreteIndividualGenerationStrategy
        : public IIndividualGenerationStrategy<IndividualType> {
      std::function<IndividualType()> generationStrategy;

    public:
      ConcreteIndividualGenerationStrategy(
          std::function<IndividualType()> generationStrategy)
          : generationStrategy(generationStrategy) {}
      IndividualType generateIndividual() override {
        return generationStrategy();
      }
    };

    this->fitnessCalculationStrategy =
        std::make_unique<ConcreteFitnessCalculationStrategy>(
            fitnessCalculationStrategy);
    this->viabilityStrategy =
        std::make_unique<ConcreteViabilityStrategy>(viabilityStrategy);
    this->individualGenerationStrategy =
        std::make_unique<ConcreteIndividualGenerationStrategy>(
            generationStrategy);
    initializePopulation();
  }

  void run(int g) {
    for (int i = 0; i < g; i++) {
      std::vector<IndividualType> elite;
      for (int j = 0; j < populationSize * elitismRate; j++) {
        elite.push_back(population[j]);
      }
      crossoverPopulation();
      mutatePopulation();
      evaluatePopulation();
      population.insert(population.end(), elite.begin(), elite.end());
      std::sort(population.begin(), population.end());

      population.erase(population.begin() + populationSize, population.end());
      // std::cout << "Generation " << i << " best fitness: " <<
      // population[0].getFitness() << std::endl;
    }
  }
  IndividualType getBestIndividual() { return population[0]; }
};
// Implementation as a Higher Order Function
// Using functional programming
#include <algorithm>
#include <functional>
#include <random>
#include <vector>

namespace FunctionalGA {

template <typename IndividualType>
class DefaultViabilityStrategy : public IViabilityStrategy<IndividualType> {
public:
  bool isViable(const IndividualType &individual) override { return true; }
};

template <typename IndividualType>
void defaultMutationFunction(IndividualType &individual) {
  std::vector<bool> genes = individual.getGenotype();
  // set up random seed
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> mut(0, 1);

  for (int j = 0; j < genes.size(); j++) {
    if (mut(gen) < 0.01) {
      genes[j] = !genes[j];
    }
  }

  individual.setGenotypeAndUpdatePhenotype(genes);
}

template <typename IndividualType>
IndividualType
defaultSelectionFunction(const std::vector<IndividualType> &population) {
  // Get a pair of distinct Individuals using tournament selection
  std::vector<IndividualType> tournament;
  // set up random seed
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, population.size() - 1);

  tournament.clear();
  for (int j = 0; j < 10; j++) {
    tournament.push_back(population[dis(gen)]);
  }
  std::sort(tournament.begin(), tournament.end());
  return tournament[0];
}

template <typename IndividualType>
std::pair<IndividualType, IndividualType>
defaultCrossoverFunction(IndividualType ind1, IndividualType ind2) {
  std::vector<bool> genes1 = ind1.getGenotype();
  std::vector<bool> genes2 = ind2.getGenotype();
  // set up random seed
  std::random_device rd;
  std::mt19937 gen(rd());

  std::uniform_int_distribution<> dis(0, genes1.size() - 1);
  int crossoverPoint = dis(gen);

  std::vector<bool> newGenes1(genes1.begin(), genes1.begin() + crossoverPoint);
  newGenes1.insert(newGenes1.end(), genes2.begin() + crossoverPoint,
                   genes2.end());
  std::vector<bool> newGenes2(genes2.begin(), genes2.begin() + crossoverPoint);
  newGenes2.insert(newGenes2.end(), genes1.begin() + crossoverPoint,
                   genes1.end());

  return std::make_pair(IndividualType(newGenes1), IndividualType(newGenes2));
}

template <typename IndividualType>
std::vector<IndividualType>
initializePopulation(std::function<IndividualType()> generationStrategy,
                     int populationSize) {
  std::vector<IndividualType> population(populationSize);
  std::generate(population.begin(), population.end(), generationStrategy);
  return population;
}

template <typename IndividualType>
std::vector<IndividualType> evaluatePopulation(
    const std::vector<IndividualType> &population,
    std::function<float(const IndividualType &)> fitnessCalculationStrategy) {
  std::vector<IndividualType> newPopulation;
  newPopulation.reserve(population.size());
  std::transform(
      population.begin(), population.end(), std::back_inserter(newPopulation),
      [fitnessCalculationStrategy](const IndividualType &individual) {
        IndividualType newIndividual = individual;
        newIndividual.setFitness(fitnessCalculationStrategy(newIndividual));
        return newIndividual;
      });
  return newPopulation;
}

template <typename IndividualType>
std::vector<IndividualType>
mutatePopulation(const std::vector<IndividualType> &population,
                 float mutationRate,
                 std::function<void(IndividualType &)> mutationFunction =
                     defaultMutationFunction<IndividualType>) {
  std::vector<IndividualType> mutatedPopulation = population;

  std::transform(mutatedPopulation.begin(), mutatedPopulation.end(),
                 mutatedPopulation.begin(),
                 [mutationRate, mutationFunction](auto &i) {
                   if (rand() / (float)RAND_MAX < mutationRate) {
                     mutationFunction(i);
                   }
                   return i;
                 });
  return mutatedPopulation;
}

template <typename IndividualType>
std::vector<IndividualType> selectParents(
    const std::vector<IndividualType> &population, int numberOfParentsToSelect,
    std::function<IndividualType(const std::vector<IndividualType> &)>
        selectionFunction = defaultSelectionFunction<IndividualType>) {
  std::vector<IndividualType> parents;
  parents.reserve(numberOfParentsToSelect);
  std::generate_n(std::back_inserter(parents), numberOfParentsToSelect,
                  [selectionFunction, &population]() {
                    return selectionFunction(population);
                  });
  return parents;
}

template <typename IndividualType>
std::vector<IndividualType> crossoverPopulation(
    const std::vector<IndividualType> &parents, int numberOfChildrenToGenerate,
    std::function<std::pair<IndividualType, IndividualType>(IndividualType,
                                                            IndividualType)>
        crossoverFunction = defaultCrossoverFunction<IndividualType>) {
  std::vector<IndividualType> newPopulation;
  newPopulation.reserve(numberOfChildrenToGenerate * 2);

  auto it = parents.begin();
  std::generate_n(
      std::back_inserter(newPopulation), numberOfChildrenToGenerate * 2, [&]() {
        if (it == parents.end()) {
          it = parents.begin();
        }
        auto firstParent = *it++;
        if (it == parents.end()) {
          it = parents.begin();
        }
        auto secondParent = *it++;
        auto children = crossoverFunction(firstParent, secondParent);
        return children.first;
      });
  return newPopulation;
}

} // namespace FunctionalGA
