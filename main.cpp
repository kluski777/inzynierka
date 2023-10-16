#include <iostream>
#include <cstdlib>
#include <cmath>

typedef long double ldouble;
const int N = 64000;
const int zElem = 4;
const int xElem = zElem*10, yElem = zElem*10;
// założenie też jest takie że wymiary xElem == yElem
// do tego liczba komórek jest zawsze parzysta

const int xElemD2 = xElem/2, yElemD2 = yElem/2, zElemD2 = zElem/2;

//////////////////////////// Dane wielkości
const ldouble T = 6.93792e15; 				// okres obrotu w sekundach
const ldouble M = 6e42;
const ldouble mMin = M/N;	 					// masa 1 cząstki w kg 
const ldouble G = 6.67e-11; 					// stała grawitacyjna w SI.
const ldouble w = 2*M_PI/T; 					// częstość kołowa rad.
const ldouble dt = T/1000;						// krok czasowy, minimum 1000 kroków na obrót
const ldouble dx = 9.258e19/zElem;			// wymiar jednej kratki

//////////////////////////// Pomocnicze żeby przyspieszyć obliczenia
const ldouble dx2 = dx*dx;
const ldouble coswt = cos(w*dt);
const ldouble sinwt = sin(w*dt);
const ldouble Gm = G*mMin;
const ldouble dt2mdx = dt*dt/(dx*mMin);
const ldouble GM = G*M;
const ldouble m2 = mMin/2;
const ldouble dt2 = 1/(2*2*dt*dt);
const ldouble potConst = 2/3*M_PI*dx2*G;
const ldouble one6th = 1.0/6;

class massPoint{ // totalna zmiana tej klasy, wszystko out.
	public:
		ldouble oldX, oldY, oldZ, x, y, z;
		massPoint(): oldX{dx*(rand()%xElem - xElemD2)}, oldY{dx*(rand()%yElem - yElemD2)}, oldZ{dx*(rand()%zElem - zElemD2)}, z{oldZ}, y{oldX*sinwt + oldY*coswt}, x{oldX*coswt - oldY*sinwt}{}
};

class dane{
	public:
		massPoint* place; 	// KONIECZNIE ZMIENIĆ NA ARRAY 3D, w którym też jest ilość
		int*** density; 		// ilość ciał w kratce
		ldouble*** pot, E;
		ldouble*** space;

		dane(); 					// warunki początkowe
		~dane();
		void potEdge(); 		// przypisanie ścian
		void potAssign();
		ldouble step();
		void printIndxOut();	// na Eulera usunąć te dwie funkcje poniżej.
		void printPotOut();
};

inline int xCoord2indx(ldouble coord){	return int(coord/dx + xElemD2); }
inline int yCoord2indx(ldouble coord){ return int(coord/dx + yElemD2); }
inline int zCoord2indx(ldouble coord){	return int(coord/dx + zElemD2); }

void verletStep(massPoint point, ldouble*** potential, int*** density, ldouble& E){ // trzeba sobie lepiej poradzić z tą energią
	ldouble X = point.x, Y = point.y, Z = point.z;

	int i = xCoord2indx(X), j = yCoord2indx(Y), k = zCoord2indx(Z);

	if( i >= 0 && i < xElem && j >= 0 && j < yElem && k >= 0 && k < zElem)
		density[i][j][k] -= 1;
	
	ldouble r3 = pow(X*X+Y*Y+Z*Z, 3.0/2);

	if(i > 0 && i < xElem-1)
		point.x = 2*X - point.oldX - (potential[i+1][j][k] - potential[i-1][j][k])*dt2mdx;
	else
		point.x = 2*X - point.oldX - GM*X/r3;
	
	if(j > 0 && j < yElem-1)
		point.y = 2*Y - point.oldY - (potential[i][j+1][k] - potential[i][j-1][k])*dt2mdx;
	else
		point.y = 2*Y - point.oldY - GM*Y/r3;

	if(k > 0 && k < zElem-1)
		point.z = 2*Z - point.oldZ - (potential[i][j][k+1] - potential[i][j][k-1])*dt2mdx;
	else
		point.z = 2*Z - point.oldZ - GM*Z/r3;

	E += dt2*(pow(point.x - point.oldX, 2) + pow(point.y - point.oldY, 2) + pow(point.z - point.oldZ, 2));	// potencjał się doda w klasie danych.

	point.oldX = X;
	point.oldY = Y;
	point.oldZ = Z;

	i = xCoord2indx(point.x);
	j = yCoord2indx(point.y);
	k = zCoord2indx(point.z);

	if( i >=0 && i < xElem && j >= 0 && j < yElem && k >= 0 && k < zElem)
		density[i][j][k] += 1;

}

int main(){
	dane symEZ;

	for(int i=0; i<1000; i++){
		ldouble E = symEZ.step();
		std::cout << i+1 << ". " << E << std::endl;
		if(std::isinf(E) || std::isnan(E)){
			printf("Iteracja %d, E = %Le\n", i, E);
			break;
		}
	}

	symEZ.printIndxOut();

	return 0;
}

ldouble dane::step(){
	E = 0;
	for(int i=0; i<N; i++)
		verletStep(place[i], pot, density, E);
	
	printf("Same prędkości: %Le\n", E);

	E *= m2;
	printf("Po verlecie E = %Le\n", E);

	this->potAssign();
	for(int i=0; i<xElem; i++){
		for(int j=0; j<yElem; j++){
			for(int k=0; k<zElem; k++)
				E += density[i][j][k]*pot[i][j][k];
		}
	}
	return E;
}

void dane::printPotOut(){
	for(int i=0; i<xElem; i++){
		for(int j=0; j<yElem; j++){
			for(int k=0; k<zElem; k++)
				printf("%Le ", pot[i][j][k]);
			printf("\n");
		}
		printf("\n");
	}
}

void dane::potAssign(){
	for(int i=1; i<xElem-1; i++){
		for(int j=1; j<yElem-1; j++){
			for(int k=1; k<zElem-1; k++)
				pot[i][j][k] = one6th*(pot[i-1][j][k]+pot[i+1][j][k]+pot[i][j-1][k]+pot[i][j+1][k]+pot[i][j][k-1]+pot[i][j][k+1]) - potConst*density[i][j][k];
		}
	}
}

void dane::printIndxOut(){
	printf("Elementy które wyrzuciło z planszy:\n");
	for(int i=0; i<N; i++){
		int xIndx = int(place[i].x/dx + xElemD2), yIndx = int(place[i].y/dx + yElemD2), zIndx = int(place[i].z/dx + zElemD2);
		if( xIndx > xElem || yIndx > yElem || zIndx > zElem || xIndx < 0 || yIndx < 0 || zIndx < 0 )
			printf("place[%d] = [%d][%d][%d]\n", i, xIndx, yIndx, zIndx);
	}
}

dane::~dane(){
	for(int x=0; x<xElem; x++){
		for(int y=0; y<yElem; y++){
			delete [] space[x][y];
			delete [] pot[x][y];
			delete [] density[x][y];
		}
		delete [] space[x];
		delete [] pot[x];
		delete [] density[x];
	}
	delete [] space;
	delete [] pot; 
	delete [] density;
	delete [] place;
}

dane::dane(){
	////////////////// alokacja pamięci na masę, potencjał oraz tablicowanie odległości 
	space = new ldouble**[xElem];
	density = new int**[xElem];
	pot = new ldouble**[xElem];
	for(int i=0; i<xElem; i++){
		density[i] = new int*[yElem];
		pot[i] = new ldouble*[yElem];
		space[i] = new ldouble*[yElem];
		for(int j=0; j<yElem; j++){
			space[i][j] = new ldouble[zElem];
			density[i][j] = new int[zElem];
			pot[i][j] = new ldouble[zElem];		
			for(int k=0; k<zElem; k++)
				space[i][j][k] = 1/sqrt(i*i*dx2 + j*j*dx2 + k*k*dx2); // to można zredukować !!!
		}
	}

	place = new massPoint[N];
	/////////////////// koniec alokacji, początek przypisywania
	for(int i=0; i<N; i++){
		place[i] = massPoint(); 
		int xIdx = int(place[i].x/dx + xElemD2), yIdx = int(place[i].y/dx + yElemD2), zIdx = int(place[i].z/dx + zElemD2);
		if( xIdx >= 0 && xIdx < xElem && yIdx >= 0 && yIdx < yElem && zIdx >= 0 && zIdx < zElem) 	// warunek żeby nie było segFaulta
			density[xIdx][yIdx][zIdx] += 1;
	}

	potEdge();
	potAssign();
}

void dane::potEdge(){ // trzeba to kiedyś ulepszyć, jak się da / może na kartce algorytm policzyć (Co ja miałem na myśli?)
	for(int i=0; i<xElem; i++){
		for(int j=0; j<yElem; j++){
			pot[i][j][0] = 0.0;
			pot[i][j][zElem-1] = 0.0;
			for(int k=0; k<zElem; k++){
				for(int wspX=0; wspX<xElem; wspX++){
					for(int wspY=0; wspY<yElem; wspY++){
						pot[i][j][0] -= i==wspX && j == wspY && k == 0 ? 0 : density[wspX][wspY][k] * space[abs(i-wspX)][abs(j-wspY)][k];
						pot[i][j][zElem-1] -= i == wspX && j == wspY && k == zElem-1 ? 0 : density[wspX][wspY][zElem - 1 - k] * space[abs(i-wspX)][abs(j-wspY)][zElem - 1 - k];
					}
				}
			}
			pot[i][j][0] *= Gm;
			pot[i][j][zElem-1] *= Gm; 
		}
	}

	// brzeg x
	for(int i=0; i<xElem; i++){ // x albo y, bez różnicy.
		for(int k=1; k<zElem-1; k++){
			pot[i][0][k] = 0;
			pot[i][yElem-1][k] = 0;	// czy zerowanie jest konieczne
			for(int wspX=0; wspX<xElem; wspX++){
				for(int wspY=0; wspY<yElem; wspY++){
					for(int wspZ=0; wspZ<zElem; wspZ++){
						pot[i][0][k] -= wspX == i && wspY == 0 && wspZ == k ? 0 : density[wspX][wspY][wspZ] * space[abs(i-wspX)][wspY][abs(wspZ-k)];
						pot[i][yElem-1][k] -= wspX == i && wspY == yElem-1 && wspZ == k ? 0 : density[wspX][wspY][wspZ] * space[abs(i-wspX)][yElem-1-wspY][abs(wspZ-k)];
					}
				}
			}
			pot[i][0][k] *= Gm;
			pot[i][yElem-1][k] *= Gm;
		}
	}
	
	// brzeg y, może by to upchać do jednej pętli
	for(int j=1; j<yElem-1; j++){
		for(int k=1; k<zElem-1; k++){
			pot[0][j][k] = 0;
			pot[xElem-1][j][k] = 0;
			for(int wspX=0; wspX<xElem; wspX++){
				for(int wspY=0; wspY<yElem; wspY++){
					for(int wspZ=0; wspZ<zElem; wspZ++){
						pot[0][j][k] -=  wspX == 0 && wspY == j && wspZ == k ? 0 : density[wspX][wspY][wspZ] * space[wspX][abs(wspY-j)][abs(wspZ-k)];
						pot[xElem-1][j][k] -= wspX == xElem-1 && wspY == j && wspZ == k ? 0 : density[wspX][wspY][wspZ] * space[xElem-1-wspX][abs(wspY-j)][abs(wspZ-k)];
					}
				}
			}
			pot[0][j][k] *= Gm;
			pot[xElem-1][j][k] *= Gm;
		}
	}
}
