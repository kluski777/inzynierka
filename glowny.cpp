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
const ldouble T = 6.93792e15; 					// okres obrotu w sekundach
const ldouble M = 6e42;
const ldouble mMin = M/N;	 					// masa 1 cząstki w kg 
const ldouble G = 6.67e-11; 					// stała grawitacyjna w SI.
const ldouble w = 2*M_PI/T; 					// częstość kołowa rad.
const ldouble dt = T/1000;						// krok czasowy, minimum 1000 kroków na obrót
const ldouble dx = 9.258e19/zElem;				// wymiar jednej kratki

//////////////////////////// Pomocnicze żeby przyspieszyć obliczenia
const ldouble dx2 = dx*dx;
const ldouble kratkaDensity = mMin/(dx2*dx); 	// jedna cząstka w kratce daje taką gęstość na kratkę
const ldouble coswt = cos(w*dt);
const ldouble sinwt = sin(w*dt);
const ldouble Gm = G*mMin;
const ldouble dt2mdx = dt*dt/(dx*mMin);
const ldouble GM = G*M;
const ldouble dt2p8 = dt*dt*8;

class massPoint{
	public:
		ldouble x, y, z, oldX, oldY, oldZ, E;
		massPoint(): oldX{dx*(rand()%xElem - xElemD2)}, oldY{dx*(rand()%yElem - yElemD2)}, oldZ{0}, z{0}, y{oldX*sinwt + oldY*coswt}, x{oldX*coswt - oldY*sinwt}, E{0} {}
};

class dane{
	public:
		massPoint* place; 	// punkty masowe, da się to jakoś lepiej rozwiązać ?
		int*** density; 	// ilość ciał w kratce
		ldouble*** pot;
		ldouble*** space;

		dane(); 			// warunki początkowe
		void potAssign(); 	// przypisanie ścian
};

void verletStep(massPoint point, ldouble*** potential, int i, int j, int k){
	ldouble X = point.x, Y = point.y, Z = point.z;

	if(i > 0 && i < xElem-1)
		point.x = 2*X - point.oldX - (potential[i+1][j][k] - potential[i-1][j][k])*dt2mdx;
	else
		point.x = 2*X - point.oldX - GM/(X*X+Y*Y+Z*Z) * cos(asin(Z/sqrt(X*X+Y*Y+Z*Z))) * X/sqrt(X*X+Y*Y);
	
	if(j > 0 && j < yElem-1)
		point.y = 2*Y - point.oldY - (potential[i][j+1][k] - potential[i][j-1][k])*dt2mdx;
	else
		point.y = 2*Y - point.oldY - GM/(X*X+Y*Y+Z*Z) * cos(asin(Z/sqrt(X*X+Y*Y+Z*Z))) * sin(acos(X/sqrt(X*X+Y*Y)));

	if(k > 0 && k < zElem-1)
		point.z = 2*Z - point.oldZ - (potential[i][j][k+1] - potential[i][j][k-1])*dt2mdx;
	else
		point.z = 2*Z - point.oldZ - GM * Z / sqrt(X*X+Y*Y+Z*Z);
	
	int xIndx = int(point.x/dx - xElemD2), yIndx = int(point.y/dx - yElemD2), zIndx = int(point.z/dx - zElemD2);
	point.E = mMin*(pot[xIndx][yIndx][zIndx] + /())

	point.oldX = X;
	point.oldY = Y;
	point.oldZ = Z;
}

int main(){
	dane symEZ;
	

	return 0;
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
		if( xIdx >= 0 && xIdx <= xElem-1 && yIdx >= 0 && yIdx < yElem-1 && zIdx >= 0 && zIdx < zElem-1) 	// warunek żeby nie było segFaulta
			density[xIdx][yIdx][zIdx] += 1;
	}

	potAssign();
	std::cout << "Gęstości:" << std::endl;
	for(int i=0; i<xElem; i++){
		for(int j=0; j<yElem; j++)
			printf("[%d][%d][2]. %d\n", i, j, density[i][j][2]);
	}
}

void dane::potAssign(){
	// dla z != 0 && zElem-1, jeszcze nie wiem czy to jest dobrze, pierwsze pętle są dla ścianek z maksymalnym z
	for(int i=0; i<xElem; i++){
		for(int j=0; j<yElem; j++){
			for(int k=1; k<zElem-1; k++){
				for(int wspX=1; wspX<xElem-1; wspX++){
					for(int wspY=1; wspY<yElem-1; wspY++){
						pot[i][j][0] -=  density[wspX][wspY][k] * space[abs(i-wspX)][abs(j-wspY)][k];
						pot[i][j][zElem-1] -= density[wspX][wspY][zElem - 1 - k] * space[abs(i-wspX)][abs(j-wspY)][zElem - 1 - k];
					}
				}
			}
			pot[i][j][0] *= Gm;
			pot[i][j][zElem-1] *= Gm; 
		}
	}

	// tj. dla pozostałych ścianek.
	for(int i=0; i<xElem; i++){ // x albo y, bez różnicy.
		for(int k=0; k<zElem; k++){
			for(int wspX=1; wspX<xElem-1; wspX++){
				for(int wspY=1; wspY<yElem-1; wspY++){
					for(int wspZ=1; wspZ<zElem-1; wspZ++){
						pot[i][0][k] -= density[wspX][wspY][wspZ] * space[abs(i-wspX)][wspY][abs(wspZ-k)];
						pot[i][yElem-1][k] -= density[wspX][wspY][wspZ] * space[abs(i-wspX)][yElem-1-wspY][abs(wspZ-k)];
						pot[0][i][k] -= density[wspX][wspY][wspZ] * space[wspX][abs(wspY-i)][abs(wspZ-k)];
						pot[xElem-1][i][k] -= density[wspX][wspY][wspZ] * space[xElem-1-wspX][abs(wspY-i)][abs(wspZ-k)];
					}
				}
			}
			pot[i][0][k] *= Gm;
			pot[i][yElem-1][k] *= Gm;
			pot[0][i][k] *= Gm;
			pot[xElem-1][i][k] *= Gm;
		}
	}
}
