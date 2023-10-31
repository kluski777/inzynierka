#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>

typedef long double ldouble;
const int N = 64000;
const int zElem = 4;
const int xElem = zElem*10, yElem = zElem*10;
	// założenie też jest takie że wymiary xElem == yElem
	// do tego liczba komórek jest zawsze parzysta

const int xElemD2 = xElem/2, yElemD2 = yElem/2, zElemD2 = zElem/2; // 20, 20, 2

//////////////////////////// Dane wielkości
const ldouble T = 6.93792e15; 				// okres obrotu w sekundach
const ldouble M = 6e42;
const ldouble mMin = M/N;	 					// masa 1 cząstki w kg 
const ldouble G = 6.67e-11; 					// stała grawitacyjna w SI.
const ldouble w = 2*M_PI/T; 					// częstość kołowa rad.
const ldouble dt = T/1000;					// krok czasowy, minimum 1000 kroków na obrót
const ldouble dx = 9.258e19/zElem;			// wymiar jednej kratki

//////////////////////////// Pomocnicze żeby przyspieszyć obliczenia
const ldouble dx2 = dx*dx;
const ldouble dx12 = .5*dx;
const ldouble coswt = cos(w*dt);
const ldouble sinwt = sin(w*dt);
const ldouble Gm = G*mMin;
const ldouble dt2mdx = dt*dt/(2*dx*mMin); // rzędu 10^-32, w mianowniku 2 powinno być
const ldouble GM = G*M;
const ldouble m2 = mMin/2;
const ldouble dt2 = dt*dt;
const ldouble dt4 = 1.0/(4*dt2);
const ldouble potConst = (2.0/3.0)*M_PI*dx2*G*mMin/pow(dx,3.0);
const ldouble one6th = 1.0/6.0;

class massPoint{
	public:
		ldouble oldX, oldY, oldZ, x, y, z;
		int liczba;
		massPoint(ldouble wspX, ldouble wspY, ldouble wspZ, int ile): oldX{wspX}, oldY{wspY}, oldZ{wspZ}, z{oldZ}, y{oldX*sinwt + oldY*coswt}, x{oldX*coswt - oldY*sinwt}, liczba{ile}{}
		massPoint(){}
};

class dane{
	public:
		massPoint*** place; 	// położenia ciał, do verleta potrzebne
		int*** density; 		// ilość ciał w kratce
		ldouble*** pot, E;
		ldouble*** space;

		dane(); 					// warunki początkowe
		~dane();
		void potEdge(); 		// przypisanie brzegów ścian potencjału
		void potAssign();		// przypisanie potencjałów dla reszty pudła
		ldouble step();
};

inline int xCoord2indx(ldouble coord){	return round(coord/dx + xElemD2); }
inline int yCoord2indx(ldouble coord){ return round(coord/dx + yElemD2); }
inline int zCoord2indx(ldouble coord){	return round(coord/dx + zElemD2); }

void verletStep(massPoint & point, ldouble*** pot, int*** density, ldouble& E){
	ldouble X = point.x;
	ldouble Y = point.y;
	ldouble Z = point.z;

	int i = xCoord2indx(X);
	int j = yCoord2indx(Y);
	int k = zCoord2indx(Z);

	if( i>=0 && i<xElem && j>=0 && j<yElem && k>=0 && k<zElem)
		density[i][j][k] -= point.liczba;

	ldouble cr3 = dt2/(pow(X*X+Y*Y+Z*Z, 3.0/2.0));

	if(i > 0 && i < xElem-1 && j>= 0 && j<yElem && k>=0 && k<zElem){
		point.x = 2*X - point.oldX - (pot[i+1][j][k] - pot[i-1][j][k])*dt2mdx;
		//printf("Wyraz 1: %Le, wyraz 2: %Le\n", 2*X - point.oldX, (pot[i+1][j][k] - pot[i-1][j][k])*dt2mdx);
	}
	else{ // 38 rzędzów wielkości różnicy między tymi dwoma sposobami. 
		point.x = 2*X - point.oldX - GM*X*cr3;
		//printf("Wyraz 1: %Le, wyraz 2: %Le\n", 2*X - point.oldX, GM*X*cr3); // 4 rzędy wielkości różnicy.
	}
	
	if(j > 0 && j < yElem-1 && i>=0 && i<xElem && k>=0 && k<zElem)
		point.y = 2*Y - point.oldY - (pot[i][j+1][k] - pot[i][j-1][k])*dt2mdx;
	else
		point.y = 2*Y - point.oldY - GM*Y*cr3;

	if(k > 0 && k < zElem-1 && i>=0 && i<xElem && j>=0 && j<yElem)
		point.z = 2*Z - point.oldZ - (pot[i][j][k+1] - pot[i][j][k-1])*dt2mdx;
	else
		point.z = 2*Z - point.oldZ - GM*Z*cr3;

	E += dt4*(pow(point.x - point.oldX, 2) + pow(point.y - point.oldY, 2) + pow(point.z - point.oldZ, 2));

	point.oldX = X;
	point.oldY = Y;
	point.oldZ = Z;

	i = xCoord2indx(point.x);
	j = yCoord2indx(point.y);
	k = zCoord2indx(point.z);
	
	if( i >= 0 && i < xElem && j >= 0 && j < yElem && k >= 0 && k < zElem)
		density[i][j][k] -= point.liczba;
}

int main(){
	dane symEZ;

	std::ofstream fileOut("toPlot.txt");
	std::ofstream pot("Potencjaly", std::ios::binary); // potem można dodać poruszające się wykresy niby.

	for(int i=0; i<10000; i++){
		printf("Energia = %Le\n", symEZ.step());
		fileOut << symEZ.place[10][10][0].x << " " << symEZ.place[10][10][0].y << " " << symEZ.place[10][10][0].z << std::endl;
	}

	for(int i=0; i<xElem; i++){
		for(int j=0; j<yElem; j++)
			pot << symEZ.pot[i][j][0] << " ";
		pot << std::endl;
	}

	return 0;
}

ldouble dane::step(){
	E = 0;
	for(int i=0; i<xElem; i++){
		for(int j=0; j<yElem; j++){
			for(int k=0; k<zElem; k++)
				verletStep(place[i][j][k], pot, density, E);
		}
	}
	
	E *= m2;

	potAssign();
	for(int i=0; i<xElem; i++){
		for(int j=0; j<yElem; j++){
			for(int k=0; k<zElem; k++)
				E += density[i][j][k]*pot[i][j][k];
		}
	}
	return E;
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
				space[i][j][k] = 1.0/(dx*sqrt(double(i*i + j*j + k*k))); // z tablicowanie odległości
		}
	}

	place = new massPoint**[xElem];
	for(int i=0; i<xElem; i++){
		place[i] = new massPoint*[yElem];
		for(int j=0; j<yElem; j++){
			place[i][j] = new massPoint[zElem];
			for(int k=0; k<zElem; k++)
				place[i][j][k] = massPoint( dx*(i-xElemD2) + dx12, dx*(j-yElemD2) + dx12, dx*(k-zElemD2) + dx12, 0);
		}
	}
	////////////// Koniec alokacji, przypisywanie.

	int xIndx, yIndx, zIndx;
	for(int i=0; i<xElem; i++){
		for(int j=0; j<yElem; j++){
			for(int k=0; k<zElem*10; k++){
				xIndx = rand()%xElem;
				yIndx = rand()%yElem;
				zIndx = rand()%zElem;
				place[xIndx][yIndx][zIndx].liczba++; // ilość cząstek w komórce
			}
		}
	}

	for(int i=0; i<xElem; i++){
		for(int j=0; j<yElem; j++){
			for(int k=0; k<zElem; k++)
				density[i][j][k] = place[i][j][k].liczba;
		}
	}

	potEdge();
}


void dane::potAssign(){ //
	for(int i=1; i<xElem-1; i++){
		for(int j=1; j<yElem-1; j++){
			for(int k=1; k<zElem-1; k++)
				pot[i][j][k] = one6th*(pot[i-1][j][k]+pot[i+1][j][k]+pot[i][j-1][k]+pot[i][j+1][k]+pot[i][j][k-1]+pot[i][j][k+1]) - potConst*density[i][j][k];
		}
	}
}


void dane::potEdge(){
	// Część z przypisywaniem ścianek do poprawy
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
	for(int i=0; i<xElem; i++){
		for(int k=1; k<zElem-1; k++){
			pot[i][0][k] = 0;
			pot[i][yElem-1][k] = 0;	// czy zerowanie jest konieczne?
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

	bool hasChanged = true;
	while(hasChanged){
		hasChanged = false;
		for(int i=1; i<xElem-1; i++){
			for(int j=1; j<yElem-1; j++){
				for(int k=1; k<zElem-1; k++){
					ldouble temp = pot[i][j][k];
					pot[i][j][k] = one6th*(pot[i-1][j][k]+pot[i+1][j][k]+pot[i][j-1][k]+pot[i][j+1][k]+pot[i][j][k-1]+pot[i][j][k+1]) - potConst*density[i][j][k];
					if( temp != pot[i][j][k] )
						hasChanged = true;
				}
			}
		}
	}
}

// TODO cząstki nie wyłażą za siatkę, ale część przyspieszenia od potencjałów jest śmiesznie niska. (potencjalnie stała do sprawdzenia)
// TODO niedokładności na poziomie ~5% przy przypisaniu wartości początkowych, ~1/6 dx
