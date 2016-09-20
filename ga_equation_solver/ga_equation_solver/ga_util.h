#ifndef GA_UTIL_H
#define GA_UTIL_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

// params of GA
const double PTYPE_MAX = 5.0;
const double PTYPE_MIN = -5.0;
const int GTYPE_LEN = 20;
const int GTYPE_MAX = 1;
const int POP = 10;
const int NUM_OF_GENERATION = 100;

// params of GA Operators
const double ELITE_RATE = 0.1;
const double MUTATE_RATE = 0.02;

// GA Farget Function
double fx(double x) {
	return pow((x), 2) - 9.0;
}

// GA Fitness Function
double gy(double y) {
	return 1.0 / (1.0 + fabs(y));
}

/*
 * Geno Type
 */
class Gtype {
public:
	Gtype() {};

	void setGtypeRamdom(int gtypeCodeLength, int gtypeCodeMax);
	int mutateGtype();
	int getGCodeLength() { return gCodeLength; }
	int getGCodeMax() { return gCodeMax; }

	Gtype operator =(Gtype& g) {
		this->gVector = g.gVector;
		this->gCodeLength = g.getGCodeLength();
		this->gCodeMax = g.getGCodeMax();
		return *this;
	}

	std::vector<int> gVector;

private:
	int gCodeLength;
	int gCodeMax;
};

void Gtype::setGtypeRamdom(int gtypeCodeLength, int gtypeCodeMax){
	gCodeLength = gtypeCodeLength;
	gCodeMax = gtypeCodeMax;
	gVector.resize(gtypeCodeLength);

	for (int i = 0; i < gtypeCodeLength; i++)
		gVector[i] = rand() % (gtypeCodeMax + 1);
}

int Gtype::mutateGtype(){
	int mutateCount = 0;
	int mutateFlag = 0;

	for (int i = 0; i < gCodeLength; i++) {
		if (((double)rand() / RAND_MAX) <= MUTATE_RATE) {
			gVector[i] = (gVector[i] == 1)? 0 : 1;
			mutateCount++;
		}
	}

	return mutateCount;
}

/*
 * Pheno Type
 */
class Ptype {
public:
	Ptype() {};

	void decodeGtype(Gtype& g);
	void setRange(double min, double max);
	void setPtype(double ptype) { p = ptype; };
	double getPtype() { return p; };

private:
	double p;
	double pMin;
	double pMax;
};

void Ptype::decodeGtype(Gtype& g) {
	double gap = pMax - pMin;
	double decoded_value = pMin;

	int i = 0;

	while (i < g.getGCodeLength()) {
		if (g.gVector[i])
			decoded_value += gap / pow(2, i+1);
		i++;
	}

	p = decoded_value;
}

void Ptype::setRange(double min, double max) {
	pMin = min;
	pMax = max;
}

/*
 * individual
 */
class GaIndividual {
public:
	GaIndividual() {};
	
	void setGaIndividual(int gtypeCodeLength, int gtypeCodeMax, double ptyeMin, double ptypeMax);
	void calcFitness(double p);

	static bool compareGaIndividualPredicate(GaIndividual a, GaIndividual b) { return (a.fitness > b.fitness); }

	double getFitness() { return fitness; }
	int getRank() { return rank; }
	int getParent1() { return parent1; }
	int getParent2() { return parent2; }
	int getCrossPoint() { return crossPoint; }
	void setRank(int gRank) { rank = gRank; }
	void setParent1(int p1) { parent1 = p1; }
	void setParent2(int p2) { parent2 = p2; }
	void setCrossPoint(int cp) { crossPoint = cp; }

	Gtype gtype; // Geno Type
	Ptype ptype; // Pheno Type

private:
	double fitness;
	int rank;  // rank after sorting
	int parent1; // index of parent 1
	int parent2; // index of parent 2
	int crossPoint; // crossover point
};

void GaIndividual::setGaIndividual(int gtypeCodeLength, int gtypeCodeMax, double ptyeMin, double ptypeMax){
	rank = 0;
	parent1 = 0;
	parent2 = 0;
	crossPoint = 0;

	gtype.setGtypeRamdom(gtypeCodeLength,gtypeCodeMax);
	ptype.setRange(ptyeMin, ptypeMax);
}

void GaIndividual::calcFitness(double p) {
	fitness = gy(fx(p));
}

/*
 * a generation
 */
class GaGenaration {
public:
	GaGenaration() {};

	void initGeneration(int gaPopulationSize, int gtypeCodeLength, int gtypeCodeMax, double ptyeMin, double ptypeMax);
	void coutResult(int num);
	int getPopulation() { return static_cast<int>(genes.size()); }

	void evaluation();
	void selection();
	void crossover();
	void mutation();

private:

	GaIndividual decideParentRoulette(std::vector<GaIndividual>& candidate);
	GaIndividual getChildrenCrossOver(GaIndividual& parent1, GaIndividual& parent2);

	std::vector<GaIndividual> genes;
	int mutateCount;  // total number of mutation
	double maxFitness;
	double minFitness;
	double avgFitness;
};

void GaGenaration::initGeneration(int gaPopulationSize, int gtypeCodeLength, int gtypeCodeMax, double ptyeMin, double ptypeMax){
	mutateCount = 0;
	maxFitness = 0.0;
	minFitness = 0.0;
	avgFitness = 0.0;
	genes.resize(gaPopulationSize);

	for (int i = 0; i < getPopulation(); i++)
		genes[i].setGaIndividual(gtypeCodeLength, gtypeCodeMax, ptyeMin, ptypeMax);
}

void GaGenaration::evaluation() {
	// calc ptype from gtype
	for (int i = 0; i < getPopulation(); i++)
		genes[i].ptype.decodeGtype(genes[i].gtype);
}

void GaGenaration::selection(){
	// calc fitness of each genes
	avgFitness = 0.0;
	for (int i = 0; i < getPopulation(); i++) {
		genes[i].calcFitness(genes[i].ptype.getPtype());
		avgFitness += genes[i].getFitness();
	}
	avgFitness /= genes.size();

	// sort genes by fitness
	std::sort(genes.begin(), genes.end(), GaIndividual::compareGaIndividualPredicate);

	// identify genes
	for (int i = 0; i < getPopulation(); i++)
		genes[i].setRank(i);

	// calc min and max of fitness
	maxFitness = genes[0].getFitness();
	minFitness = genes[getPopulation() - 1].getFitness();
}

void GaGenaration::coutResult(int num) {

	std::cout << "----------------------------------------------------" << std::endl;
	std::cout << "###  parent   xsite  gtype                    ptype   fitness" << std::endl;
	int count = 0;
	for (auto gene : genes) {
		std::cout << std::setw(3) << count << " (" << std::setw(3) << gene.getParent1() << "," << std::setw(3) << gene.getParent2() << ")  " << std::setw(3) << gene.getCrossPoint() << "   [";
		for (int i = 0; i < gene.gtype.getGCodeLength(); i++)
			std::cout << gene.gtype.gVector[i];
		std::cout << "] " << std::setw(6) << gene.ptype.getPtype() << "    " << std::setw(6) << gene.getFitness() << std::endl;
		count++;
	}
	std::cout << std::endl << "Mutation: " << mutateCount << std::endl
		<< num << ", " << maxFitness << ", " << avgFitness << ", " << minFitness
		<< ", " << genes[0].ptype.getPtype() << ", [";
	for (int i = 0; i < genes[0].gtype.getGCodeLength(); i++)
		std::cout << genes[0].gtype.gVector[i];
	std::cout << "]" << std::endl;
}

void GaGenaration::crossover(){
	// Reincarnation of Elite
	for (int i = 0; i < getPopulation()*ELITE_RATE; i++) {
		genes[i].setParent1(0);
		genes[i].setParent2(0);
		genes[i].setCrossPoint(0);
	}

	// decide parent
	std::vector<GaIndividual> parent1;
	std::vector<GaIndividual> parent2;
	GaIndividual p2Candidate;
	for (int i = 0; i < getPopulation()*(1.0 - ELITE_RATE); i++) {
		parent1.push_back(decideParentRoulette(genes));
		p2Candidate = decideParentRoulette(genes);
		while (p2Candidate.getRank() == parent1[i].getRank()) {
			p2Candidate = decideParentRoulette(genes);
		}
		parent2.push_back(p2Candidate);
	}

	//std::cout << std::endl << parent1.size() << ", " << parent1.size() << std::endl;

	// cross over parent
	std::vector<GaIndividual> children;
	for (int i = 0; i < getPopulation()*(1.0 - ELITE_RATE); i++)
		children.push_back(getChildrenCrossOver(parent1[i], parent2[i]));

	// modify generation
	int numOfChildren = static_cast<int>(children.size());
	for (int i = 0; i < numOfChildren; i++)
		genes[getPopulation() - children.size() + i] = children[i];
}

void GaGenaration::mutation(){
	mutateCount = 0;
	for (int i = 0; i < getPopulation()*(1.0 - ELITE_RATE); i++)
		mutateCount += genes[getPopulation()*ELITE_RATE +i].gtype.mutateGtype();
}

GaIndividual GaGenaration::decideParentRoulette(std::vector<GaIndividual>& candidate){
	double r = (double)rand() / RAND_MAX; // 0~1
	std::vector<double> fitnessCategory;
	double totalFitness;
	int numOfCandidate = static_cast<int>(candidate.size());

	// calc fitness category and total fitness
	fitnessCategory.push_back(0.0);
	for (int i = 0; i < numOfCandidate; i++)
		fitnessCategory.push_back(fitnessCategory[i] + genes[i].getFitness());
	totalFitness = fitnessCategory[numOfCandidate];

	// decide parent
	for (int i = 0; i < numOfCandidate; i++)
		if ((r >= fitnessCategory[i] / totalFitness) && (r < fitnessCategory[i + 1] / totalFitness))
			return candidate[i];

	return GaIndividual();
}

GaIndividual GaGenaration::getChildrenCrossOver(GaIndividual& parent1, GaIndividual& parent2){
	GaIndividual child;
	child = parent1;

	//set parent data
	child.setParent1(parent1.getRank());
	child.setParent2(parent2.getRank());

	//decide intersection ramdomely
	int intersection;
	intersection = (int)(parent1.gtype.getGCodeLength() * rand() / RAND_MAX);
	child.setCrossPoint(intersection);

	//calc new gtype
	for (int i = 0; i < child.gtype.getGCodeLength(); i++)
		if (i < intersection)
			child.gtype.gVector[i] = parent2.gtype.gVector[i];

	return child;
}

#endif // GA_UTIL_H