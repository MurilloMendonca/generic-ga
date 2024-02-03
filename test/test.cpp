#include "../include/GA.hpp"
#include "../include/Helper.hpp"
#include "../include/Individual.hpp"
#include <chrono>
#include <functional>
#include <iostream>
#include <random>
using namespace std;

bool testGrayConversion() {
  float a = 0.5463;
  std::vector<bool> gray = Helper::float2gray(a, 32);

  cout << "decimal: " << a << endl;
  cout << "gray: ";
  for (auto i : gray) {
    cout << i;
  }

  float b = Helper::gray2float(gray);

  cout << "\ndecimal: " << b << endl;
  return abs(a - b) < 0.0001;
}
class IntGenotypePhenotypeStrategy : public GenotypePhenotypeStrategy<int> {
public:
  int genotype2phenotype(const std::vector<bool> &genotype) override {

    return Helper::gray2decimal(genotype);
  }

  std::vector<bool> phenotype2genotype(const int &phenotype) override {
    // Implement conversion logic here
    return Helper::decimal2gray(phenotype, 10);
  }
};

class FloatGenotypePhenotypeStrategy : public GenotypePhenotypeStrategy<float> {
public:
  float genotype2phenotype(const std::vector<bool> &genotype) override {
    std::vector<bool> decimalPart(genotype.begin(), genotype.begin() + 10);
    std::vector<bool> floatPart(genotype.begin() + 10, genotype.end());
    int dec = Helper::gray2decimal(decimalPart);
    float floaty = Helper::gray2float(floatPart);
    return dec + floaty;
  }

  std::vector<bool> phenotype2genotype(const float &phenotype) override {
    // Implement conversion logic here
    auto phenotypeInt = Helper::decimal2gray(phenotype, 10);
    auto phenotypeFloat = Helper::float2gray(phenotype - (int)phenotype, 30);
    phenotypeInt.insert(phenotypeInt.end(), phenotypeFloat.begin(),
                        phenotypeFloat.end());
    return phenotypeInt;
  }
};

class ConcreteFitnessStrategy
    : public IFitnessCalculationStrategy<Individual<float>> {
public:
  float calculateFitness(const Individual<float> &individual) override {
    // Implement fitness calculation logic here
    const float x = individual.getPhenotype();
    return sin(x) * sin(2 * x + 5);
  }
};

class ConcreteViabilityStrategy : public IViabilityStrategy<Individual<float>> {
public:
  bool isViable(const Individual<float> &individual) override {
    // Implement viability check logic here
    const float x = individual.getPhenotype();
    return x >= 0 && x <= 10;
  }
};

class ConcreteIndividualGenerationStrategy
    : public IIndividualGenerationStrategy<Individual<float>> {
public:
  Individual<float> generateIndividual() override {
    // Implement individual generation logic here
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 10);
    return Individual<float>(dis(gen));
  }
};

bool testOOPWithInterfaces() {
  FloatGenotypePhenotypeStrategy *strategy =
      new FloatGenotypePhenotypeStrategy();
  Individual<float>::setStrategy(strategy);

  auto fitness = std::make_unique<ConcreteFitnessStrategy>();
  auto isViable = std::make_unique<ConcreteViabilityStrategy>();
  auto generateValidIndividual =
      std::make_unique<ConcreteIndividualGenerationStrategy>();
  int generations = 100;
  int populationSize = 100;
  GA<Individual<float>> ga(std::move(fitness), std::move(isViable),
                           std::move(generateValidIndividual), populationSize,
                           0.1, 0.1);
  Individual<float> bestIndividual = ga.getBestIndividual();
  for (int i = 0; i < generations; i++) {
    ga.run(1);
    auto genBestIndividual = ga.getBestIndividual();
    bestIndividual = std::max(bestIndividual, genBestIndividual,
                              [](const auto &lhs, const auto &rhs) {
                                return lhs.getFitness() < rhs.getFitness();
                              });
  }

  std::cout << "Best individual: " << bestIndividual.getPhenotype() << " - "
            << bestIndividual.getFitness() << std::endl;

  return bestIndividual.getFitness() > 0.98;
}

bool testOOPWithLambdas() {
  auto individualGenerator = []() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 10);
    return Individual<float>(dis(gen));
  };
  auto fitnessCalc = [](const Individual<float> &individual) {
    const float x = individual.getPhenotype();
    return sin(x) * sin(2 * x + 5);
  };
  auto isViableFunc = [](const Individual<float> &individual) {
    const float x = individual.getPhenotype();
    return x >= 0 && x <= 10;
  };

  GA<Individual<float>> ga(fitnessCalc, isViableFunc, individualGenerator, 100,
                           0.1, 0.1);
  auto bestIndividual = ga.getBestIndividual();
  for (int i = 0; i < 100; i++) {
    ga.run(1);
    auto genBestIndividual = ga.getBestIndividual();
    bestIndividual = std::max(bestIndividual, genBestIndividual,
                              [](const auto &lhs, const auto &rhs) {
                                return lhs.getFitness() < rhs.getFitness();
                              });
  }

  std::cout << "Best individual: " << bestIndividual.getPhenotype() << " - "
            << bestIndividual.getFitness() << std::endl;
  return bestIndividual.getFitness() > 0.98;
}

bool testFunctional() {
  auto individualGenerator = []() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 10);
    return Individual<float>(dis(gen));
  };
  auto fitnessCalc = [](const Individual<float> &individual) {
    const float x = individual.getPhenotype();
    return sin(x) * sin(2 * x + 5);
  };
  auto isViableFunc = [](const Individual<float> &individual) {
    const float x = individual.getPhenotype();
    return x >= 0 && x <= 10;
  };

  auto population = FunctionalGA::initializePopulation<Individual<float>>(
      individualGenerator, 100);

  auto bestIndividual =
      *std::max_element(population.begin(), population.end(),
                        [](const auto &lhs, const auto &rhs) {
                          return lhs.getFitness() < rhs.getFitness();
                        });
  for (int i = 1; i < 100; i++) {
    population = FunctionalGA::evaluatePopulation<Individual<float>>(
        population, fitnessCalc);
    auto genBestIndividual =
        *std::max_element(population.begin(), population.end(),
                          [](const auto &lhs, const auto &rhs) {
                            return lhs.getFitness() < rhs.getFitness();
                          });
    bestIndividual =
        genBestIndividual.getFitness() > bestIndividual.getFitness()
            ? genBestIndividual
            : bestIndividual;
    population = FunctionalGA::mutatePopulation<Individual<float>>(
        FunctionalGA::crossoverPopulation<Individual<float>>(
            FunctionalGA::selectParents<Individual<float>>(population, 100 / 3),
            100),
        1.1);
  }

  std::cout << "Best individual: " << bestIndividual.getPhenotype() << " - "
            << bestIndividual.getFitness() << std::endl;
  return bestIndividual.getFitness() > 0.98;
}
int main() {
  bool result = true;
  std::cout << "============ Gray Conversion Test ============" << std::endl;
  auto start = std::chrono::high_resolution_clock::now();
  result = testGrayConversion();
  auto end = std::chrono::high_resolution_clock::now();
  if (result) {
    std::cout << "Test passed" << std::endl;
  } else {
    std::cout << "Test failed" << std::endl;
  }
  std::cout << "Time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                     start)
                   .count()
            << "ms" << std::endl;
  std::cout << "============ OOP with Interfaces Test ============"
            << std::endl;
  start = std::chrono::high_resolution_clock::now();
  result = testOOPWithInterfaces();
  end = std::chrono::high_resolution_clock::now();
  if (result) {
    std::cout << "Test passed" << std::endl;
  } else {
    std::cout << "Test failed" << std::endl;
  }
  std::cout << "Time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                     start)
                   .count()
            << "ms" << std::endl;
  std::cout << "============ OOP with Lambdas Test ============" << std::endl;
  start = std::chrono::high_resolution_clock::now();
  result = testOOPWithLambdas();
  end = std::chrono::high_resolution_clock::now();
  if (result) {
    std::cout << "Test passed" << std::endl;
  } else {
    std::cout << "Test failed" << std::endl;
  }
  std::cout << "Time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                     start)
                   .count()
            << "ms" << std::endl;
  std::cout << "============ Functional Test ============" << std::endl;
  start = std::chrono::high_resolution_clock::now();
  testFunctional();
  end = std::chrono::high_resolution_clock::now();
  if (result) {
    std::cout << "Test passed" << std::endl;
  } else {
    std::cout << "Test failed" << std::endl;
  }
  std::cout << "Time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                     start)
                   .count()
            << "ms" << std::endl;
  return 0;
}
