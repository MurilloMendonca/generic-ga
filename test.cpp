#include <iostream>
#include <string>
#include <functional>
#include <random>
#include "Individual.hpp"
#include "GA.hpp"
#include "Helper.hpp"


using namespace std;

void testGrayConversion(){
	float a = 0.5463;
	string gray = Helper::float2gray(a, 25);

	cout << "decimal: " << a << endl;
	cout << "gray: " << gray << endl;

	float b = Helper::gray2float(gray);

	cout << "decimal: " << b << endl;
}
class IntGenotypePhenotypeStrategy : public GenotypePhenotypeStrategy<int> {
public:
    int genotype2phenotype(const std::string& genotype) override {
        
        return Helper::gray2decimal(genotype);
    }

    std::string phenotype2genotype(const int& phenotype) override {
        // Implement conversion logic here
        return Helper::decimal2gray(phenotype, 10);
    }
};

class FloatGenotypePhenotypeStrategy : public GenotypePhenotypeStrategy<float> {
public:
    float genotype2phenotype(const std::string& genotype) override {
        
        int decimalPart = Helper::gray2decimal(genotype.substr(0, 10));
		float floatPart = Helper::gray2float(genotype.substr(10, 30));
		return decimalPart + floatPart;
    }

    std::string phenotype2genotype(const float& phenotype) override {
        // Implement conversion logic here
        return Helper::decimal2gray(phenotype, 10) + Helper::float2gray(phenotype - (int)phenotype, 20);
    }
};


int main()
{
	//testGrayConversion();
	FloatGenotypePhenotypeStrategy* strategy = new FloatGenotypePhenotypeStrategy();
	Individual<float>::setStrategy(strategy);

	std::function<float(Individual<float>)> fitness = [](Individual<float> i) {
		float x = i.getPhenotype();
		return sin(x)*sin(2*x+5);
	};

	std::function<bool(Individual<float>)> isViable = [](Individual<float> i) {
		float x = i.getPhenotype();
		return x>=0 && x<=10;
	};

	std::function<Individual<float>()> generateValidIndividual = []() {
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<> dis(0, 10);
		return Individual<float>(dis(gen));
	};

	GA<Individual<float>> ga(generateValidIndividual, fitness, isViable, 10, 0.1, 0.1);

	for(int i = 0; i < 100; i++)
	{
		ga.run(1);
		cout << ga.getBestIndividual().getPhenotype() << endl;
	}


	
	return 0;
}