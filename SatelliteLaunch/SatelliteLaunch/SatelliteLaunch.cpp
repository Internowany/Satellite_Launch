/*
Project no. 5 - Satellite Launch
created by Maciej Bekas & Sebastian Kalinowski

Given a simple (two-dimensional) model of a planetary system, prepare a program that given an
origin planet, calculates the time, direction vector and speed needed to launch a satellite, so
that it flies close to a given destination planet.
Assume that we are only interested in the gravitational interactions of the planets with the
satellite and that all planets have the same mass.
Input:
• Initial (time 𝑡 = 0) positions of the planets (we may assume circular orbits)
  (𝑟1 ,𝜙1),… , (𝑟𝑁,𝜙𝑁) in polar coordinates and their orbital periods 𝑝1,… ,𝑝𝑁,
• the index of the origin and destination planets,
• the distance within which the satellite needs pass next to the destination planet,
• the maximum allowed flight time.
Output: time, angle, and satellite launch speed.
*/

#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;

/*
USER INPUTS
*/
const int iterations = 500;
int sourcePlanet = 4;
int destinationPlanet = 1;
double distancePlanetSatellite = 0.5;	// million km
const int missionTime = 600;			// days
const int numberOfPlanets = 4;			// we MUST assume number of planets

//Parameters for mutation
double par1 = 0.5;						// for speed
double par2 = 2.0;						// for angle

// Initial positions of the planets (radius, angle, period) are defined in main() function

/*
CONSTANTS
*/
double maxAngle = 2 * M_PI;
double maxVelocity = 3.6288;			// million km / day
double planetR = 0.012;					// million km
double GM = 2.97545 * pow(10, -3);		// (million km)^3 / day^2

/*
CLASSES AND DATA STRUCTURES
*/
struct Tablica							// single dimensional array - Array[2]
{
	double first, second;
};
class Planet
{
public:
	Planet(double r, double angle, double period);

	int getId();
	double getR();
	double getAngle();
	double getPeriod();
	Planet* getNext();

	void setId(int id);
	void setR(double r);
	void setAngle(double angle);
	void setPeriod(double period);
	void setNext(Planet* next);

	double calculateAngle(double time);
	friend ostream& operator<<(ostream& os, const Planet& planet);
private:
	int id;
	double r, angle, period;
	Planet* next;
};
class System
{
public:
	System();
	~System();

	Planet* getHead();

	int addPlanet(double radius, double angle, double period);
	void sortPlanets();
	void swap(Planet* a, Planet* b);
	Planet* searchPlanet(int id);
	double calculatePosition(int id, double deltaT);
	friend ostream& operator<<(ostream& os, const System& system);
private:
	Planet* head;
};
class Satellite
{
public:
	Satellite(Planet* planet, double t0, double velocity, double angle);

	double getR();
	double getAngle();
	double getVr();
	double getT0();
	Tablica getXy();
	Tablica getVxvy();

	void setXy(Tablica xy);
	void setVxvy(Tablica vxvy);

	Tablica calculatePosition(double time);
private:
	double r, angle, vr, t0;
	Tablica xy, vxvy;
};
class Simulation
{
public:
	Simulation(System* system, Satellite* satellite);

	void updateSatellite(int time, double interval);
	Tablica simulateFlight(int missionTime, Planet* destination);
private:
	System* system;
	Satellite* satellite;
};
double tab[missionTime+1][14];
double tab2[missionTime+1][14];

/*
FUNCTIONS
*/
double fRand(double fMin, double fMax);
Tablica convToCartesian(double r, double angle);
double calculateDistance(double x1, double y1, double x2, double y2);
double generateGaussianNoise(double mu, double sigma);
double mutate(double x, double sigma);
double HillClimbing(System* sys);
void saveToCsv(double bestTime, const char* csvFile);

// Planet methods
Planet::Planet(double r, double angle, double period)
{
	this->id = 0;
	this->r = r;
	this->angle = angle;
	this->period = period;
	this->next = NULL;
}
int Planet::getId()
{
	return id;
}
double Planet::getR()
{
	return r;
}
double Planet::getAngle()
{
	return angle;
}
double Planet::getPeriod()
{
	return period;
}
Planet* Planet::getNext()
{
	return next;
}
void Planet::setId(int id)
{
	this->id = id;
}
void Planet::setR(double r)
{
	this->r = r;
}
void Planet::setAngle(double angle)
{
	this->angle = angle;
}
void Planet::setPeriod(double period)
{
	this->period = period;
}
void Planet::setNext(Planet* next)
{
	this->next = next;
}
double Planet::calculateAngle(double time)
{
	double the = this->getAngle() + time / this->getPeriod() * maxAngle;
	return fmod(the, maxAngle);
}
ostream& operator<<(ostream& os, const Planet& planet)
{
	os << "Id: " << planet.id << "; Radius: " << planet.r << "; Delta: " << planet.angle << "; Period: " << planet.period << endl;
	return os;
}

// System methods
System::System()
{
	head = NULL;
	cout << "Planet system:\n\n";
}
System::~System()
{
	cout << "\nPlanets destroyed!\n-------------------------------------------------------------\n";
	delete head;
}
Planet* System::getHead()
{
	return head;
}
int System::addPlanet(double r, double angle, double period)
{
	Planet* current = head;
	Planet* previous;

	while (current != NULL)
	{
		if (r == current->getR())  //check if there is a planet with the same id
		{
			cout << "Error: Planet with this radius already exists:" << r << endl << endl;
			return 1;
		}
		current = current->getNext();
	}
	current = head;

	Planet* new_planet;
	new_planet = new Planet(r, angle, period);

	if (head != NULL)
	{
		while (current != NULL)
		{
			if (current->getNext() == NULL)    //to add if it is last element
			{
				current->setNext(new_planet);
				new_planet->setNext(NULL);
				return 0;
			}
			previous = current;
			current = current->getNext();
		}
	}
	else if (head == NULL)  //to add first element
	{
		new_planet->setNext(head);
		head = new_planet;
		return 0;
	}
	return 1;
}
void System::sortPlanets()
{
	int swapped, i = 1;
	Planet* current;
	Planet* next = NULL;

	do
	{
		swapped = 0;
		current = head;

		while (current->getNext() != next)
		{
			if (current->getR() > current->getNext()->getR())
			{
				swap(current, current->getNext());
				swapped = 1;
			}
			current = current->getNext();
		}
		next = current;

	} while (swapped);

	current = head;
	while (current != NULL)
	{
		current->setId(i);
		i++;
		current = current->getNext();
	}
}
void System::swap(Planet* a, Planet* b)
{
	double temp1 = a->getR();
	double temp2 = a->getAngle();
	double temp3 = a->getPeriod();

	a->setR(b->getR());
	a->setAngle(b->getAngle());
	a->setPeriod(b->getPeriod());

	b->setR(temp1);
	b->setAngle(temp2);
	b->setPeriod(temp3);
}
Planet* System::searchPlanet(int id)
{
	Planet* curr = head;
	while (curr != NULL)
	{
		if (curr->getId() == id)
		{
			return curr;
		}
		curr = curr->getNext();
	}
	return NULL;
}
double System::calculatePosition(int id, double deltaT)
{
	Planet* a = searchPlanet(id);
	double position = fmod((a->getAngle() + 2 * (deltaT / a->getPeriod()) * M_PI), (2 * M_PI));
	cout << position << endl;
	return position;
}
ostream& operator<<(ostream& os, const System& system)
{
	Planet* current = system.head;

	while (current != NULL)
	{
		cout << *current;
		current = current->getNext();
	}
	cout << endl;
	return os;
}

// Satellite methods
Satellite::Satellite(Planet* planet, double t0, double velocity, double angle)
{
	this->r = planet->getR() + planetR;
	this->angle = planet->getAngle();
	this->vr = velocity;
	this->t0 = t0;
	this->xy = convToCartesian(this->r, this->angle);
	this->vxvy = convToCartesian(this->vr, angle);
}
double Satellite::getR()
{
	return r;
}
double Satellite::getAngle()
{
	return angle;
}
double Satellite::getVr()
{
	return vr;
}
double Satellite::getT0()
{
	return t0;
}
Tablica Satellite::getXy()
{
	return xy;
}
Tablica Satellite::getVxvy()
{
	return vxvy;
}
void Satellite::setXy(Tablica xy1)
{
	this->xy.first = xy1.first;
	this->xy.second = xy1.second;
}
void Satellite::setVxvy(Tablica vxvy1)
{
	this->vxvy.first = vxvy1.first;
	this->vxvy.second = vxvy1.second;
}
Tablica Satellite::calculatePosition(double time)
{
	Tablica position;
	position.first = this->xy.first + this->vxvy.first * time;
	position.second = this->xy.second + this->vxvy.second * time;
	return position;
}

// Simulation methods
Simulation::Simulation(System* sys, Satellite* sat)
{
	this->system = sys;
	this->satellite = sat;
}
void Simulation::updateSatellite(int time, double interval)
{
	Planet* curr = system->getHead();
	Tablica axay;
	axay.first = 0;
	axay.second = 0;
	int b = 1;
	Tablica sxsy = this->satellite->getXy();
	tab[time + 1][1] = sxsy.first;
	tab[time + 1][2] = sxsy.second;

	while (curr != NULL)
	{
		Tablica pxpy = convToCartesian(curr->getR(), curr->calculateAngle(time));
		double distance = calculateDistance(this->satellite->getXy().first, this->satellite->getXy().second, pxpy.first, pxpy.second);
		//cout << "distance: " << distance << endl;
		double cos = (this->satellite->getXy().first - pxpy.first) / distance;
		double sin = (this->satellite->getXy().second - pxpy.second) / distance;
		double a = GM / pow(distance, 2);
		axay.first += a * cos;
		axay.second += a * sin;
		//cout << a <<endl;
		//cout << axay.first << " " << axay.second << endl;
		curr = curr->getNext();
		tab[time + 1][2 * b + 1] = pxpy.first;
		tab[time + 1][2 * b + 2] = pxpy.second;
		b++;
	}
	Tablica vsxvsy;
	//cout << satellite->getVxvy().first << " " << satellite->getVxvy().second << endl;
	vsxvsy.first = this->satellite->getVxvy().first - axay.first * interval;
	vsxvsy.second = this->satellite->getVxvy().second - axay.second * interval;
	//cout << axay.first << " " << axay.second << " " <<vsxvsy.first << " " << vsxvsy.second << endl;
	//cout << this->satellite->getVxvy().first << " " << this->satellite->getVxvy().second << endl;
	this->satellite->setVxvy(vsxvsy);
	sxsy = this->satellite->calculatePosition(interval);

	this->satellite->setXy(sxsy);
}
Tablica Simulation::simulateFlight(int missionTime, Planet* destination)
{
	double bestDist = 9999999999999;
	int time;
	for (time = 0; time < missionTime; time++)
	{
		if (time >= this->satellite->getT0())
		{
			this->updateSatellite(time, 1);
		}
		Tablica txty = convToCartesian(destination->getR(), destination->calculateAngle(time));
		double distance = calculateDistance(this->satellite->getXy().first, this->satellite->getXy().second, txty.first, txty.second);
		if (bestDist > distance)
		{
			bestDist = distance;
		}
		if (bestDist < distancePlanetSatellite)
		{
			//cout << bestDist;
			break;
		}
	}
	Tablica wynik;
	wynik.first = time;
	wynik.second = bestDist;
	return wynik;
}

// Misc functions
double fRand(double fMin, double fMax)
{
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}
Tablica convToCartesian(double r, double angle)
{
	Tablica cartCord;
	cartCord.first = r * cos(angle);
	cartCord.second = r * sin(angle);
	return cartCord;
}
double calculateDistance(double x1, double y1, double x2, double y2)
{
	return sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2));
}
double generateGaussianNoise(double mu, double sigma)	// we used Box-Muller transform for generation of Gaussian Noise used in mutation https://en.wikipedia.org/wiki/Box–Muller_transform
{
	static const double epsilon = numeric_limits<double>::min();

	thread_local double z1;
	thread_local bool generate;
	generate = !generate;

	if (!generate)
		return z1 * sigma + mu;

	double u1, u2;
	do
	{
		u1 = rand() * (1.0 / RAND_MAX);
		u2 = rand() * (1.0 / RAND_MAX);
	} while (u1 <= epsilon);

	double z0;
	z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
	z1 = sqrt(-2.0 * log(u1)) * sin(2.0 * M_PI * u2);

	return z0 * sigma + mu;
}
double mutate(double x, double sigma)
{
	return x + generateGaussianNoise(0, sigma);
}
double HillClimbing(System* sys)
{
	/* Hill climbing algorithm v.2
	x - speed/angle; t - speedMutated/angleMutated; f(x) - bestDist; f(t) - bestDistMutated; par1/par2 - mutation parameters

	Stop conditions:
	- max iterations
	- satellite reached a destination
	*/
	double speed = maxVelocity / (par1 * iterations);
	double angle = maxAngle / (par2 * iterations);
	double bestDist = 9999999999999;
	double bestTime = 0;

	for (int i = 0; i < iterations; i++)
	{
		Planet* target = sys->searchPlanet(destinationPlanet);

		double speedMutated = mutate(speed, maxVelocity / (par1 * iterations));		// Here we can change parameters of mutation of speed and angle, because this two values 
		double angleMutated = mutate(angle, maxAngle / (par2 * iterations));		// have different units, so it is wise to mutate those values with different parameters.
		if (angleMutated < 0)
		{
			angleMutated = 2 * M_PI + angleMutated;
		}
		if (angleMutated > 2 * M_PI)
		{
			angleMutated = angleMutated - 2 * M_PI;
		}
		if (speedMutated > 0)
		{
			Satellite* satellite2 = new Satellite(sys->searchPlanet(sourcePlanet), 0, speed, angle);	// In each iteration simulate for x and t, then compare f(x) and f(t)
			Simulation* simulation2 = new Simulation(sys, satellite2);
			Tablica bestDist = simulation2->simulateFlight(missionTime, target);

			Satellite* satellite = new Satellite(sys->searchPlanet(sourcePlanet), 0, speedMutated, angleMutated);
			Simulation* simulation = new Simulation(sys, satellite);
			Tablica bestDistMutated = simulation->simulateFlight(missionTime, target);


			if (bestDistMutated.second < bestDist.second)		// if f(t) < f(x)
			{
				speed = speedMutated;							// x <- t
				angle = angleMutated;
				for (int i = 0; i <= missionTime; i++)
				{
					for (int j = 1; j <= 14; j++)
					{
						tab2[i][j] = tab[i][j];
					}
				}
				bestTime = bestDistMutated.first;
			}
		}
	}
	cout << "Best speed: " << speed << " million km / day" << endl;
	cout << "Best angle: " << angle << " radians" << endl;
	cout << "Best time: " << bestTime << " days" << endl;
	return bestTime;
}
void saveToCsv(double bestTime, const char* csvFile)
{
	int i, j;
	FILE* fp;

	if ((fp = fopen(csvFile, "w")) == NULL)
	{
		printf("ERROR!\n");
		exit(1);
	}
	for (j = 1; j <= bestTime + 1; j++)
	{
		for (i = 1; i <= 2 * numberOfPlanets + 2; i = i + 2)
		{
			fprintf(fp, "%f,", tab2[j][i]);
		}
		fprintf(fp, "\n");
		for (i = 2; i <= 2 * numberOfPlanets + 2; i = i + 2)
		{
			fprintf(fp, "%f,", tab2[j][i]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

int main()
{
	System* system1 = new System();		// system1
	system1->addPlanet(58, 1, 88);
	system1->addPlanet(108, 2, 224);
	system1->addPlanet(190, 3.5, 365);
	system1->addPlanet(228, 4, 686);
	system1->sortPlanets();
	cout << *system1;
	saveToCsv(HillClimbing(system1), "data.csv");
	delete system1;

	System* system2 = new System();		// system2
	system2->addPlanet(48, 1, 88);
	system2->addPlanet(178, 2, 224);
	system2->addPlanet(190, 0, 365);
	system2->addPlanet(278, 4, 686);
	system2->sortPlanets();
	cout << *system2;
	saveToCsv(HillClimbing(system2), "data2.csv");
	delete system2;

	return 0;
}