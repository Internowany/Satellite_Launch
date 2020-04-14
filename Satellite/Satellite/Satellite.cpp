#define _USE_MATH_DEFINES

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;

double bleeeeee = 23;
double maxAngle = 2 * M_PI;
double maxVelocity = 1.44;
double missionTime = 500;
double maxDelay = missionTime;
double planetR = 0.006;
int iterations = 50;
double GM = 0.0029752;
double a = 23;
double AniajestFajna = 69;
dsadassad
dsadassadsda

sdadas


struct Tablica
{
	double first, second;
};
class Planet
{
public:
	Planet(double r, double angle);

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
private:
	int id;
	double r, angle, period;
	Planet* next;
};
class System
{
public:
	System(int numOfPlanets, double factor);
	~System();

	Planet* getHead();

	int addPlanet(double radius, double angle);
	void printPlanets();
	void sortPlanets();
	void swap(Planet* a, Planet* b);
	Planet* searchPlanet(int id);
	double calculatePosition(int id, double deltaT);
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

	void updateSatellite(double time, double interval);
	double simulateFlight(int missionTime, Planet* destination);
private:
	System* system;
	Satellite* satellite;
};

double fRand(double fMin, double fMax);
Tablica convToCartesian(double r, double angle);
double calculateDistance(Satellite* satellite, Planet* planet, double time);
double calculateDistancePolar(double r1, double angle1, double r2, double angle2);
double calculateDistanceCart(double x1, double y1, double x2, double y2);
double randomGaussian(double mu, double sigma);
double genAngleStep();
double genSpeedStep();
double genDelayStep();
double mutate(double x, double sigma);

//Planet
Planet::Planet(double r, double angle)
{
	this->id = 0;
	this->r = r;
	this->angle = angle;
	double c = 0.06;
	this->period = sqrt(c * pow(r, 3));
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

//System
System::System(int numOfPlanets, double factor)
{
	head = NULL;
	cout << "Planet system:\n\n";

	srand(time(NULL));
	double angle;
	double r = 50;

	double first = fRand(0, maxAngle); //pierwsza wartosc pominieta!

	for (int i = 0; i < numOfPlanets; i++)
	{
		double a_random_double = fRand(0, maxAngle);
		angle = round(a_random_double * 1000.0) / 1000.0;
		this->addPlanet(r, angle);
		r = round(r * factor);
	}
	this->sortPlanets();
}
System::~System()
{
	cout << "Planets destroyed!\n";
	delete head;
}
Planet* System::getHead()
{
	return head;
}
int System::addPlanet(double radius, double angle)
{
	Planet* current = head;
	Planet* previous;

	while (current != NULL)
	{
		if (radius == current->getR())  //check if there is planet with the same id
		{
			cout << "Error: Planet with this radius already exists:" << radius << endl << endl;
			return 1;
		}
		current = current->getNext();
	}
	current = head;

	Planet* new_planet;
	new_planet = new Planet(radius, angle);

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
void System::printPlanets()
{
	Planet* current = head;

	while (current != NULL)
	{
		cout << "Id: " << current->getId() << "; Radius: " << current->getR() << "; Delta: " << current->getAngle() << "; Period: " << current->getPeriod() << endl;
		current = current->getNext();
	}
	cout << endl;
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

//Satellite
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

//Simulation
Simulation::Simulation(System* sys, Satellite* sat)
{
	this->system = sys;
	this->satellite = sat;
}
void Simulation::updateSatellite(double time, double interval)
{
	Planet* curr = system->getHead();
	Tablica axay;
	axay.first = 0;
	axay.second = 0;

	while (curr != NULL)
	{
		Tablica pxpy = convToCartesian(curr->getR(), curr->calculateAngle(time));
		double distance = calculateDistanceCart(this->satellite->getXy().first, this->satellite->getXy().second, pxpy.first, pxpy.second);
		double cos = (this->satellite->getXy().first - pxpy.first) / distance;
		double sin = (this->satellite->getXy().second - pxpy.second) / distance;
		double a = GM / pow(distance, 2);
		axay.first += a * cos;
		axay.second += a * sin;

		curr = curr->getNext();
	}
	Tablica vsxvsy;
	vsxvsy.first = this->satellite->getVxvy().first - (axay.first) * pow(interval, 2) / 2;
	vsxvsy.second = this->satellite->getVxvy().second - (axay.second) * pow(interval, 2) / 2;

	this->satellite->setVxvy(vsxvsy);
	Tablica sxsy = this->satellite->calculatePosition(interval);
	this->satellite->setXy(sxsy);
}
double Simulation::simulateFlight(int missionTime, Planet* destination)
{
	double bestDist = 9999999999999999;
	for (int time = 1; time < missionTime; time++)
	{
		if (time >= this->satellite->getT0())
		{
			this->updateSatellite(time, 1);
		}
		Tablica txty = convToCartesian(destination->getR(), destination->calculateAngle(time));
		double distance = calculateDistanceCart(this->satellite->getXy().first, this->satellite->getXy().second, txty.first, txty.second);
		if (bestDist > distance)
		{
			bestDist = distance;
		}
	}
	return bestDist;
}

//Misc
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
double calculateDistance(Satellite* satellite, Planet* planet, double time)
{
	double planet_angle = planet->getAngle() + 2 * M_PI / planet->getPeriod() * time;
	double a = (satellite->getR() * cos(satellite->getAngle()) + satellite->getVr() * cos(satellite->getAngle()) - planet->getR() * cos(planet_angle));
	double b = (satellite->getR() * sin(satellite->getAngle()) + satellite->getVr() * sin(satellite->getAngle()) - planet->getR() * sin(planet_angle));
	return sqrt(pow(a, 2) + pow(b, 2));
}
double calculateDistancePolar(double r1, double angle1, double r2, double angle2)
{
	return sqrt(pow(r1, 2) + pow(r2, 2) - 2 * r1 * r2 * cos(angle2 - angle1));
}
double calculateDistanceCart(double x1, double y1, double x2, double y2)
{
	return sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2));
}
double randomGaussian(double mu, double sigma)
{
	double first = fRand(0.0, 1.0); //pierwsza wartosc pominieta!
	double r1 = fRand(0.0, 1.0);
	double r2 = fRand(0.0, 1.0);
	double z1 = sqrt(-2 * log(r1)) * sin(2 * M_PI * r2);
	double z2 = sqrt(-2 * log(r1)) * cos(2 * M_PI * r2);
	double x1 = mu + z1 * sigma;
	double x2 = mu + z1 * sigma;
	return x1;
}
double genAngleStep()
{
	return maxAngle / iterations;
}
double genSpeedStep()
{
	return maxVelocity / iterations;
}
double genDelayStep()
{
	return maxDelay / iterations;
}
double mutate(double x, double sigma)
{
	return x + randomGaussian(0, sigma);
}

int main()
{
	System* system1 = new System(5, 1.618);
	system1->printPlanets();

	for (int i = 1; i < iterations; i++)
	{
		;;
	}
	return 0;
}