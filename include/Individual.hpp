#include <vector>
template <typename T> class GenotypePhenotypeStrategy {
public:
  virtual T genotype2phenotype(const std::vector<bool> &genotype) = 0;
  virtual std::vector<bool> phenotype2genotype(const T &phenotype) = 0;
  virtual ~GenotypePhenotypeStrategy() {}
};

template <typename T> class Individual {
private:
  T phenotype;
  std::vector<bool> genotype;
  double fitness;
  static GenotypePhenotypeStrategy<T> *strategy;

public:
  Individual() { this->fitness = 0; }
  Individual(const T &nn) {
    this->phenotype = nn;
    this->genotype = strategy->phenotype2genotype(nn);
    this->fitness = 0;
  }
  Individual(const std::vector<bool> &genotype) {
    this->setGenotypeAndUpdatePhenotype(genotype);
    this->fitness = 0;
  }
  void setFitness(float fitness) { this->fitness = fitness; }
  float getFitness() const { return this->fitness; }
  void setGenotypeAndUpdatePhenotype(const std::vector<bool> &genotype) {
    this->genotype = genotype;
    // Logic to update phenotype
    this->phenotype = strategy->genotype2phenotype(genotype);
  }
  std::vector<bool> getGenotype() {
    if (genotype.size() == 0) {
      this->genotype = strategy->phenotype2genotype(this->phenotype);
    }
    return genotype;
  }
  T getPhenotype() const { return this->phenotype; }
  void setPhenotype(T phenotype) { this->phenotype = phenotype; }

  bool operator<(const Individual &other) const {
    return this->fitness > other.fitness;
  }
  bool operator==(const Individual &other) const {
    return this->fitness == other.fitness;
  }

  static void setStrategy(GenotypePhenotypeStrategy<T> *newStrategy) {
    strategy = newStrategy;
  }

  static GenotypePhenotypeStrategy<T> *getStrategy() { return strategy; }
};

template <typename T>
GenotypePhenotypeStrategy<T> *Individual<T>::strategy = nullptr;
