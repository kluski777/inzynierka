#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include "stale_i_funkcje.hpp"

class massPoint{
	public:
		ldouble oldX, oldY, oldZ, x, y, z;
		int liczba;
		massPoint(ldouble wspX, ldouble wspY, ldouble wspZ, int ile): oldX{wspX}, oldY{wspY}, oldZ{wspZ}, z{oldZ}, y{oldX*sinwt + oldY*coswt}, x{oldX*coswt - oldY*sinwt}, liczba{ile}{} // tego sie chyba nie da lepiej napisać
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
		void potAss1stIter();
		void savePot();
		ldouble step(int);
};

void verletStep(massPoint & point, ldouble*** pot, int*** density, ldouble& E){ // to i potAssign przerobić
	ldouble X = point.x;
	ldouble Y = point.y;
	ldouble Z = point.z;

	int i = xCoord2indx(X);
	int j = yCoord2indx(Y);
	int k = zCoord2indx(Z);

	if( i>=0 && i<xElem && j>=0 && j<yElem && k>=0 && k<zElem)
		density[i][j][k] -= point.liczba;

	ldouble cr3 = dt2/(pow(X*X+Y*Y+Z*Z, 3.0/2.0));

	if(i > 0 && i < xElem-1 && j>= 0 && j<yElem && k>=0 && k<zElem)
		point.x = 2*X - point.oldX - (pot[i+1][j][k] - pot[i-1][j][k])*dt2dx;
	else
		point.x = 2*X - point.oldX - GM*X*cr3;
	
	if(j > 0 && j < yElem-1 && i>=0 && i<xElem && k>=0 && k<zElem)
		point.y = 2*Y - point.oldY - (pot[i][j+1][k] - pot[i][j-1][k])*dt2dx;
	else
		point.y = 2*Y - point.oldY - GM*Y*cr3;

	if(k > 0 && k < zElem-1 && i>=0 && i<xElem && j>=0 && j<yElem)
		point.z = 2*Z - point.oldZ - (pot[i][j][k+1] - pot[i][j][k-1])*dt2dx;
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
		density[i][j][k] += point.liczba;
}









int main(){
	dane symEZ;
	
	std::ofstream particle1("czastka1", std::ios::binary);
	std::ofstream particle2("czastka2", std::ios::binary);
	std::ofstream particle3("czastka3", std::ios::binary);
	std::ofstream particle4("czastka4", std::ios::binary);
	std::ofstream ener("Energia", std::ios::binary);

	for(int i=0; i<1000; i++){
		printf("Energia: %Le\n", symEZ.step(i));
		particle1 << symEZ.place[0][0][0].x << " " << symEZ.place[0][0][0].y << " " << symEZ.place[0][0][0].z <<  std::endl;
		particle2 << symEZ.place[10][10][0].x << " " << symEZ.place[10][10][0].y << " " << symEZ.place[10][10][0].z << std::endl;
		particle3 << symEZ.place[20][20][1].x << " " << symEZ.place[20][20][1].y << " " << symEZ.place[20][20][1].z << std::endl;
		particle4 << symEZ.place[xEdge][yEdge][zEdge].x << " " << symEZ.place[xEdge][yEdge][zEdge].y << " " << symEZ.place[xEdge][yEdge][zEdge].z << std::endl;
	}

	particle1.close();
	particle2.close();
	particle3.close();
	particle4.close();

	symEZ.savePot();

	return 0;
}









ldouble dane::step(int i){
	E = 0;
	for(int i=0; i<xElem; i++){
		for(int j=0; j<yElem; j++){
			for(int k=0; k<zElem; k++)
				verletStep(place[i][j][k], pot, density, E);
		}
	}
	
	E *= m2;

	if( i%10 == 0 )
		potEdge();
	potAssign();
	for(int k=0; k<zElem; k++){
		for(int i=0; i<xElem; i++){
			E += density[i][i][k]*pot[i][i][k];
			for(int j=i+1; j<yElem-i; j++)
				E += density[i][j][k]*pot[i][j][k] + density[j][i][k]*pot[j][i][k];
		}
	}
	return E;
}

void dane::savePot(){
	for(int k=0; k<zElem; k++){
		std::ofstream potFile("Potencjal" + std::to_string(k+1), std::ios::binary);
		for(int i=0; i<xElem; i++){
			for(int j=0; j<yElem; j++)
				potFile << pot[i][j][k] << " ";
			potFile << std::endl;
		}
	}
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
				space[i][j][k] = 1.0/(dx*sqrt(double(i*i + j*j + k*k)));
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
	for(int i=0; i<N; i++){
		xIndx = rand()%xElem;
		yIndx = rand()%yElem;
		zIndx = rand()%zElem;
		place[xIndx][yIndx][zIndx].liczba++; // ilość cząstek w komórce
	}

	for(int i=0; i<xElem; i++){
		for(int j=0; j<yElem; j++){
			for(int k=0; k<zElem; k++)
				density[i][j][k] = place[i][j][k].liczba;
		}
	}

	potEdge();			// Przypisanie do brzegów (powinno sie szybciej liczyć).
	potAss1stIter();	// Reszta tablicy liczona do zbieżności.
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

void dane::potAssign(){ // zmienić algorytm
	for(int repeat=0; repeat<3; repeat++){
		for(int k=1; k<zEdge; k++){
	  		for(int i=1; i<xEdge; i++){
				pot[i][i][k] = one6th*(pot[i-1][i][k]+pot[i+1][i][k]+pot[i][i-1][k]+pot[i][i+1][k]+pot[i][i][k+1]+pot[i][i][k-1] - potConst*density[i][i][k]);
				for(int j=i+1; j<yEdge; j++){			  
					pot[i][j][k] = one6th*(pot[i-1][j][k]+pot[i+1][j][k]+pot[i][j-1][k]+pot[i][j+1][k]+pot[i][j][k-1]+pot[i][j][k+1]) - potConst*density[i][j][k];
					pot[j][i][k] = one6th*(pot[j-1][i][k]+pot[j+1][i][k]+pot[j][i+1][k]+pot[j][i-1][k]+pot[j][i][k+1]+pot[j][i][k-1] - potConst*density[j][i][k]);
				}
			}
		}
	}
}

void dane::potEdge(){
	for(int i=0; i<xElem; i++){
		for(int j=0; j<yElem; j++){
			pot[i][j][0] = 0.0;
			pot[i][j][zEdge] = 0.0;
			for(int wspZ=0; wspZ<zElem; wspZ++){
				for(int wspX=0; wspX<xElem; wspX++){
					for(int wspY=0; wspY<yElem; wspY++){
						pot[i][j][0] -= i==wspX && j == wspY && wspZ == 0 ? 0 : density[wspX][wspY][wspZ] * space[abs(i-wspX)][abs(j-wspY)][wspZ];
						pot[i][j][zEdge] -= i == wspX && j == wspY && wspZ == zEdge ? 0 : density[wspX][wspY][zEdge - wspZ] * space[abs(i-wspX)][abs(j-wspY)][zEdge - wspZ];
					}
				}
			}
			pot[i][j][0] *= Gm;
			pot[i][j][zEdge] *= Gm; 
		}
	}

	for(int k=1; k<zEdge; k++){
		pot[0][0][k] = 0.0;
		pot[0][yEdge][k] = 0.0;
		pot[xEdge][0][k] = 0.0;
		pot[xEdge][yEdge][k] = 0.0;
		
		for(int wspX=0; wspX<xElem; wspX++){
			for(int wspY=0; wspY<yElem; wspY++){
				for(int wspZ=0; wspZ<k; wspZ++){ 
					pot[0][0][k] -= density[wspX][wspY][wspZ] * space[wspX][wspY][k-wspZ];
					pot[0][yEdge][k] -= density[wspX][wspY][wspZ] * space[wspX][yEdge-wspY][k-wspZ];
					pot[xEdge][0][k] -= density[wspX][wspY][wspZ] * space[xEdge-wspX][wspY][k-wspZ];
					pot[xEdge][yEdge][k] -= density[wspX][wspY][wspZ] * space[xEdge-wspX][yEdge-wspY][k-wspZ];
				}
				for(int wspZ=k+1; wspZ<zEdge; wspZ++){
					pot[0][0][k] -= density[wspX][wspY][wspZ] * space[wspX][wspY][wspZ-k];
					pot[0][yEdge][k] -= density[wspX][wspY][wspZ] * space[wspX][yEdge-wspY][wspZ-k];
					pot[xEdge][0][k] -= density[wspX][wspY][wspZ] * space[xEdge-wspX][wspY][wspZ-k];
					pot[xEdge][yEdge][k] -= density[wspX][wspY][wspZ] * space[xEdge-wspX][yEdge-wspY][wspZ-k];
				}
			}
		}
		pot[0][0][k] *= Gm;
		pot[xEdge][0][k] *= Gm;
		pot[0][yEdge][k] *= Gm;
		pot[xEdge][yEdge][k] *= Gm;
	}

	for(int k=1; k<zEdge; k++){
		for(int i=1; i<xEdge; i++){
			pot[0][i][k] = 0.0;
			pot[xEdge][i][k] = 0.0;
			pot[i][0][k] = 0.0;
			pot[i][yEdge][k] = 0.0;
			for(int wspX=0; wspX<xElem; wspX++){
				for(int wspY=0; wspY<yElem; wspY++){
					for(int wspZ=0; wspZ<k; wspZ++){
						pot[0][i][k] -= density[wspX][wspY][wspZ] * space[wspX][abs(i-wspY)][k-wspZ];
						pot[xEdge][i][k] -= density[wspX][wspY][wspZ] * space[xEdge-wspX][abs(i-wspY)][k-wspZ];
						pot[i][0][k] -= density[wspX][wspY][wspZ] * space[abs(wspX-i)][wspY][k-wspZ];
						pot[i][yEdge][k] -= density[wspX][wspY][wspZ] * space[abs(wspX-i)][yEdge-wspY][k-wspZ];
					}
					for(int wspZ=k+1; wspZ<zElem; wspZ++){
						pot[0][i][k] -= density[wspX][wspY][wspZ] * space[wspX][abs(i-wspY)][wspZ-k];
						pot[xEdge][i][k] -= density[wspX][wspY][wspZ] * space[xEdge-wspX][abs(i-wspY)][wspZ-k];
						pot[i][0][k] -= density[wspX][wspY][wspZ] * space[abs(wspX-i)][wspY][wspZ-k];
						pot[i][yEdge][k] -= density[wspX][wspY][wspZ] * space[abs(wspX-i)][wspY][wspZ-k];
					}
				}
			}
			for(int wspB=0; wspB<xElem; wspB++){
				for(int wspA=0; wspA<i; wspA++){
					pot[i][0][k] -= density[wspA][wspB][k] * space[i-wspA][wspB][0];
					pot[i][yEdge][k] -= density[wspA][wspB][k] * space[i-wspA][yEdge - wspB][0];
					pot[0][i][k] -= density[wspB][wspA][k] * space[wspB][i-wspA][k];
					pot[xEdge][i][k] -= density[wspB][wspA][k] * space[xEdge-wspB][i-wspA][0];
				}
				for(int wspA=i+1; wspA<yElem; wspA++){
					pot[i][0][k] -= density[wspA][wspB][k] * space[wspA-i][wspB][0];
					pot[i][yEdge][k] -= density[wspA][wspB][k] * space[wspA-i][yEdge-wspB][0];
					pot[0][i][k] -= density[wspB][wspA][k] * space[wspB][wspA-i][0];
					pot[xEdge][i][k] -= density[wspB][wspA][k] * space[xEdge-wspB][wspA-i][0];
				}
			}
			for(int wsp=0; wsp<xEdge; wsp++){
				pot[i][0][k] -= density[i][wsp+1][k] * space[0][wsp+1][0];
				pot[i][yEdge][k] -= density[i][wsp][k] * space[0][yEdge-wsp][0];
				pot[0][i][k] -= density[wsp+1][i][k] * space[wsp+1][0][0];
				pot[xEdge][i][k] -= density[wsp][i][k] * space[xEdge-wsp][0][0];
			}
			pot[0][i][k] *= Gm;
			pot[xEdge][i][k] *= Gm;
			pot[i][0][k] *= Gm;
			pot[i][yEdge][k] *= Gm;
		}
	}

}

void dane::potAss1stIter(){
	bool hasChanged = true;
	while(hasChanged){
		hasChanged = false;
		for(int i=1; i<xEdge; i++){
			for(int j=1; j<yEdge; j++){
				for(int k=1; k<zEdge; k++){
					ldouble temp = pot[i][j][k];
					pot[i][j][k] = one6th*(pot[i-1][j][k]+pot[i+1][j][k]+pot[i][j-1][k]+pot[i][j+1][k]+pot[i][j][k-1]+pot[i][j][k+1]) - potConst*density[i][j][k];
					if( fabs(temp-pot[i][j][k]) > 1e-10 )
						hasChanged = true;
				}
			}
		}
	}
}

// TODO jakiś refactoring czy coś w ten deseń żeby wszystko było bardziej czytelne
// TODO niedokładności na poziomie ~5% przy przypisaniu wartości początkowych, ~1/6 dx
