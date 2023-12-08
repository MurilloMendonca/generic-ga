#pragma once

#include <vector>
#include <algorithm>
#include <functional>
#include <random>
#include <thread>
template <typename IndividualType>
class GA{
	protected:
		virtual IndividualType generateValidIndividual() = 0;
		virtual float calculateFitness(IndividualType)=0;
		virtual bool isViable(IndividualType)=0;
		std::vector<IndividualType> population;
		float mutationRate;
		float elitismRate;
		int populationSize;

		void initializePopulation(){
			for(int i = 0; i < populationSize; i++)
			{
				population.push_back(this->generateValidIndividual());
			}
			evaluatePopulation();
			std::sort(population.begin(), population.end());
		}
		void evaluatePopulation(){
			for(auto& individual : population)
			{
				if(!this->isViable(individual))
				{
					individual.setFitness(0);
					continue;
				}
				individual.setFitness(calculateFitness(individual));
			}
		}

		void evaluatePopulationInParallel(){
			std::vector<std::thread> threads;
			for(auto& individual : population)
			{
				threads.push_back(std::thread([](IndividualType& individual, GA* ga){
					if(!ga->isViable(individual))
					{
						individual.setFitness(0);
						return;
					}
					
					individual.setFitness(ga->calculateFitness(individual));
				}, std::ref(individual), this));
			}
			for(auto& thread : threads)
			{
				thread.join();
			}
		}

		

		void mutatePopulation(){
			for(auto& individual : population)
			{
				mutation(individual);
			}
		}
		void crossoverPopulation(){
			std::vector<IndividualType> newPopulation;
			while(newPopulation.size() < populationSize)
			{
				std::pair<IndividualType, IndividualType> parents = selection();
				std::pair<IndividualType, IndividualType> children = crossover(parents.first, parents.second);
				newPopulation.push_back(children.first);
				newPopulation.push_back(children.second);
			}
			population = newPopulation;
		}

		std::pair<IndividualType, IndividualType> selection(){
			//Get a pair of distinct Individuals using tournament selection
			std::pair<IndividualType, IndividualType> parents;
			std::vector<IndividualType> tournament;
			//set up random seed
			std::random_device rd;
			std::mt19937 gen(rd());
			std::uniform_int_distribution<> dis(0, populationSize - 1);

			tournament.clear();
			for(int j = 0; j < 10; j++)
			{
				tournament.push_back(population[dis(gen)]);
			}
			std::sort(tournament.begin(), tournament.end());
			parents.first = tournament[0];
			parents.second = tournament[1];
			int i = 2;
			while(parents.first == parents.second && i < 10){
				parents.second = tournament[i++];
			}

			return parents;
		}

		std::pair<IndividualType, IndividualType> crossover(IndividualType ind1, IndividualType ind2, int numCrossoverPoints=1){
			std::string genes1 = ind1.getGenotype();
			std::string genes2 = ind2.getGenotype();
			//set up random seed
			std::random_device rd;
			std::mt19937 gen(rd());
			std::string newGenes1, newGenes2;
			int min = genes1.size() < genes2.size() ? genes1.size() : genes2.size();

			// Generate unique crossover points and sort them
			std::vector<int> crossoverPoints;
			std::uniform_int_distribution<> dis(1, min - 1); // avoiding 0 and last index
			while (crossoverPoints.size() < numCrossoverPoints) {
				int point = dis(gen);
				if (std::find(crossoverPoints.begin(), crossoverPoints.end(), point) == crossoverPoints.end()) {
					crossoverPoints.push_back(point);
				}
			}
			std::sort(crossoverPoints.begin(), crossoverPoints.end());

			// Perform n-point crossover
			int lastPoint = 0;
			bool switchParent = false;
			for (int point : crossoverPoints) {
				if (switchParent) {
					newGenes1 += genes2.substr(lastPoint, point - lastPoint);
					newGenes2 += genes1.substr(lastPoint, point - lastPoint);
				} else {
					newGenes1 += genes1.substr(lastPoint, point - lastPoint);
					newGenes2 += genes2.substr(lastPoint, point - lastPoint);
				}
				switchParent = !switchParent;
				lastPoint = point;
			}
			// Add remaining genes
			newGenes1 += switchParent ? genes2.substr(lastPoint) : genes1.substr(lastPoint);
			newGenes2 += switchParent ? genes1.substr(lastPoint) : genes2.substr(lastPoint);

			

			return std::make_pair(IndividualType(newGenes1), IndividualType(newGenes2));
		}
		void mutation(IndividualType& i){
			std::string genes = i.getGenotype();
			//set up random seed
			std::random_device rd;
			std::mt19937 gen(rd());
			std::uniform_real_distribution<> mut(0, 1);

			for(int j = 0; j < genes.size(); j++)
			{
				if(mut(gen) < mutationRate)
				{
					genes[j] = genes[j] == '0' ? '1' : '0';
				}
			}

			i.setGenotypeAndUpdatePhenotype(genes);
			
		}

	public:
		
		GA( int populationSize, float mutationRate, float elitismRate)
			{
				this->populationSize = populationSize;
				this->mutationRate = mutationRate;
				this->elitismRate = elitismRate;

				//initializePopulation();
			}

		GA()
			{
				this->populationSize = 100;
				this->mutationRate = 0.01;
				this->elitismRate = 0.1;

				//initializePopulation();
			}

		void run(int g){
			for(int i = 0; i < g; i++)
			{
				std::vector<IndividualType> elite;
				for(int j = 0; j < populationSize * elitismRate; j++)
				{
					elite.push_back(population[j]);
				}
				crossoverPopulation();
				mutatePopulation();
				evaluatePopulation();
				population.insert(population.end(), elite.begin(), elite.end());
				std::sort(population.begin(), population.end());
				
				population.erase(population.begin() + populationSize, population.end());
				//std::cout << "Generation " << i << " best fitness: " << population[0].getFitness() << std::endl;
			}
		}
		IndividualType getBestIndividual(){
			return population[0];
		}

};