#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include "stale_i_funkcje.hpp"

class massPoint{
	public:
		ldouble oldX, oldY, oldZ, x, y, z;
		int liczba;
		massPoint(ldouble wspX, ldouble wspY, ldouble wspZ, int ile): oldX{wspX}, oldY{wspY}, oldZ{wspZ}, z{oldZ}, liczba{ile}{}
		massPoint(){}
};

class dane{
	public:
		massPoint*** place; 	// położenia ciał, do verleta potrzebne
		int*** density; 		// ilość ciał w kratce
		ldouble*** pot, T, V;
		ldouble*** space;

		dane(); 					// warunki początkowe
		~dane();
		void potEdge(); 		// przypisanie brzegów ścian potencjału
		void potAssign();		// przypisanie potencjałów dla reszty pudła
		void potAss1stIter();
		void savePot();
		void step(int, ldouble*);
		void getAllInfo();
};

void verletStep(massPoint & point, ldouble*** pot, int*** density, ldouble& T, ldouble& V){
	ldouble X = point.x;
	ldouble Y = point.y;
	ldouble Z = point.z;

	int i = xCoord2indx(X);
	int j = yCoord2indx(Y);
	int k = zCoord2indx(Z);

	if( i>=0 && i<xElem && j>=0 && j<yElem && k>=0 && k<zElem)
		density[i][j][k] -= point.liczba;	// zmiana liczby cząstek w kratce.

	ldouble r = sqrt(X*X+Y*Y+Z*Z);			// odległość od środka układu.
	ldouble cr3 = dt2/(r*r*r);

	if(i > 0 && i < xElem-1 && j>= 0 && j<yElem && k>=0 && k<zElem)
		point.x = 2*X - point.oldX - (pot[i+1][j][k] - pot[i-1][j][k])*dt2dx;
	else
		point.x = 2*X - point.oldX - GM*X*cr3;
	
	if(j > 0 && j < yElem-1 && i>=0 && i<xElem && k>=0 && k<zElem)
		point.y = 2*Y - point.oldY - (pot[i][j+1][k] - pot[i][j-1][k])*dt2dx;
	else
		point.y = 2*Y - point.oldY - GM*Y*cr3;

// Jak odkomentowac linijki z dołu to nie działa, energia nie zachowana, cząstki wylatują.
//	if(k > 0 && k < zElem-1 && i>=0 && i<xElem && j>=0 && j<yElem)
//		point.z = 2*Z - point.oldZ - (pot[i][j][k+1] - pot[i][j][k-1])*dt2dx;
//	else
//		point.z = 2*Z - point.oldZ - GM*Z*cr3;

	T += point.liczba*dt4*(pow(point.x - point.oldX, 2) + pow(point.y - point.oldY, 2));

	point.oldX = X;
	point.oldY = Y;
//	point.oldZ = Z;

	i = xCoord2indx(point.x);
	j = yCoord2indx(point.y);
	k = zCoord2indx(point.z);
	
	if( i >= 0 && i < xElem && j >= 0 && j < yElem && k >= 0 && k < zElem)
		density[i][j][k] += point.liczba;
	else
		V -= GM*point.liczba/r;
}




int main(){
	dane symEZ;
//	symEZ.getAllInfo();

	std::ofstream prtcles("droga", std::ios::binary);
	std::ofstream ener("Energie" , std::ios::binary);

	ldouble* Es = new ldouble[2]; // energie, Es[0] - kinet. Es[1] - potencja.

	for(int step=0; step<1000; step++){
		symEZ.step(step, Es);
		
		ener << step << " " << Es[0] << " " << Es[1] << " " << Es[0] + Es[1] << std::endl;
		std::cout << step << " T = " << Es[0] << ", V = " << Es[1] << ", E = " << Es[0] + Es[1] << std::endl;
		for(int i=0; i<3; i++){
			for(int j=0; j<3; j++){
				for(int k=0; k<3; k++)
				prtcles << symEZ.place[i][j][k].x << " " << symEZ.place[i][j][k].y << " " << symEZ.place[i][j][k].z << std::endl;
			}
		}
		prtcles << "\t";
	}
	
	int suma=0, sumX=0, sumY=0, sumZ=0; // sumy cząstek które wywaliło z planszy

	for(int i=0; i<xElemD2; i++){
		for(int j=0; j<yElemD2; j++){
			for(int k=0; k<zElemD2; k++){
				int x = xCoord2indx(symEZ.place[i][j][k].x);
				int y = yCoord2indx(symEZ.place[i][j][k].y);
				int z = zCoord2indx(symEZ.place[i][j][k].z);
				
				suma += symEZ.place[i][j][k].liczba;

				if( x < 0 || x > xEdge)
					sumX += symEZ.place[i][j][k].liczba;
				if( y < 0 || y > yEdge)
					sumY += symEZ.place[i][j][k].liczba;
				if( z < 0 || z > zEdge)
					sumZ += symEZ.place[i][j][k].liczba;
			}
		}
	}

	printf("sum = %d, sumX = %d, sumY = %d, sumZ = %d\n", suma, sumX, sumY, sumZ);

	ener.close();
	prtcles.close();

	//symEZ.savePot();

	return 0;
}









void dane::step(int i, ldouble* Es){
	Es[0] = 0.0, Es[1] = 0.0;		// 0 to kinetyczna energia, 1 to potencjalna
	for(int k=0; k<zElemD2; k++){ 
		for(int i=0; i<xElemD2; i++){
			for(int j=0; j<yElemD2; j++)
				verletStep(place[i][j][k], pot, density, Es[0], Es[1]);
		}
	}

	Es[0] *= m2;

	if( i%3 == 0 ) // co ile liczyć ściany?
		potEdge();
	potAssign();
	for(int k=0; k<zElem; k++){
		for(int i=0; i<xElem; i++){
			Es[1] += density[i][i][k] * pot[i][i][k];
			for(int j=i+1; j<yElem; j++)
				Es[1] += density[i][j][k] * pot[i][j][k] + density[j][i][k] * pot[j][i][k];
		}
	}

	Es[1] *= mMin;
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
				space[i][j][k] = 1.0/(dx*sqrt(double(i*i + j*j + k*k))); // ztablicowane odległości.
		}
	}

	place = new massPoint**[xElem];
	for(int i=0; i<xElem; i++){
		place[i] = new massPoint*[yElem];
		for(int j=0; j<yElem; j++){
			place[i][j] = new massPoint[zElem];
			for(int k=0; k<zElem; k++)
				place[i][j][k] = massPoint( dx*(i-xElemD2)+dx12, dx*(j-yElemD2)+dx12, dx*(k-zElemD2)+dx12, 0); // duże niedokładności tu są.
		}
	}
	////////////// Koniec alokacji, przypisywanie.


	int xIndx, yIndx, zIndx;
	for(int i=0; i<N; i++){
		xIndx = rand()%xElem;
		yIndx = rand()%yElem;
		zIndx = rand()%zElem;
		place[xIndx][yIndx][zIndx].liczba++;
		density[xIndx][yIndx][zIndx]++;
	}

	// dobrze byłoby wrzucić nowe omegi tutaj
	// algorytm przerobić żeby tylko 1/8 kwadratu brał
	float max = indxDstce2(xElem, yElem, zElem); // dystans w kwadracie ineksów daje.
	float rIndx = max;
	float min = 0; // najmniejszy promień który jeszcze nie był rozważany
	while( min != max ){
		rIndx = max;
		for(int i=0; i<xElem; i++){
			for(int j=0; j<yElem; j++){
				for(int k=0; k<zElem; k++){
					float temp = indxDstce2(i, j, k);
					if( temp < rIndx && temp - min > 0.1)
						rIndx = temp;
				}
			}
		}
		int ileCzastek = 0;
		min = rIndx;
		for(int i=0; i<xElem; i++){
			for(int j=0; j<yElem; j++){
				for(int k=0; k<zElem; k++){
					float temp = indxDstce2(i, j, k);
					if(temp < min)
						ileCzastek += place[i][j][k].liczba;
				}
			}
		}
		
		ldouble w = sqrt(Gm*ileCzastek/pow(rIndx*dx2, 1.5));
		ldouble coswt = cos(dt*w);
		ldouble sinwt = sin(dt*w);
		for(int i=0; i<xElem; i++){
			for(int j=0; j<yElem; j++){
				for(int k=0; k<zElem; k++){
					if( fabs(indxDstce2(i, j, k) - rIndx) < 0.1 ){ // czy jest możliwość że się nie zakręcą?
						place[i][j][k].x = place[i][j][k].oldX*coswt - place[i][j][k].oldY*sinwt;
						place[i][j][k].y = place[i][j][k].oldX*sinwt + place[i][j][k].oldY*coswt;
					}
				}
			}
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
	for(int repeat=0; repeat<8; repeat++){
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

void dane::potEdge(){ // sprawdzić czy ściany sie zgadzają.
	for(int i=0; i<xElem; i++){
		for(int j=0; j<yElem; j++){
			pot[i][j][0] = 0.0;
			pot[i][j][zEdge] = 0.0;
			for(int wspZ=0; wspZ<zElem; wspZ++){
				for(int wspX=0; wspX<xElem; wspX++){
					for(int wspY=0; wspY<j; wspY++){
						pot[i][j][0] -= density[wspX][wspY][wspZ] * space[abs(i-wspX)][j-wspY][wspZ];
						pot[i][j][zEdge] -= density[wspX][wspY][wspZ] * space[abs(i-wspX)][j-wspY][zEdge - wspZ];
					}
					for(int wspY=j+1; wspY<yElem; wspY++){
						pot[i][j][0] -= density[wspX][wspY][wspZ] * space[abs(i-wspX)][wspY-j][wspZ];
						pot[i][j][zEdge] -= density[wspX][wspY][wspZ] * space[abs(i-wspX)][wspY-j][zEdge - wspZ];
					}
				}
				for(int wspX=0; wspX<i; wspX++){ // wspY = j
					pot[i][j][0] -= density[wspX][j][wspZ] * space[i-wspX][0][wspZ];
					pot[i][j][zEdge] -= density[wspX][j][wspZ] * space[i-wspX][0][zEdge - wspZ];
				}
				for(int wspX=i+1; wspX<xElem; wspX++){
					pot[i][j][0] -= density[wspX][j][wspZ] * space[wspX-i][0][wspZ];
					pot[i][j][zEdge] -= density[wspX][j][wspZ] * space[wspX-i][0][zEdge - wspZ];
				}
			}
			for(int wspZ=1; wspZ<zElem; wspZ++){ // dla wspY = j i wspX = i 
				pot[i][j][0] -= density[i][j][wspZ] * space[0][0][wspZ];
				pot[i][j][zEdge] -= density[i][j][wspZ] * space[0][0][zEdge - wspZ + 1];
			}
			pot[i][j][0] *= Gm;
			pot[i][j][zEdge] *= Gm;
		}
	}

	for(int k=1; k<zEdge; k++){ // nie jestem pewny czy to dziala
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


void dane::getAllInfo(){
	std::ofstream brzegDol("DolnyPot", 		std::ios::binary);
	std::ofstream brzegGora("GornyPot", 	std::ios::binary);
	std::ofstream brzegPrawy("PrawyPot", 	std::ios::binary);
	std::ofstream brzegLewy("LewyPot", 		std::ios::binary);
	std::ofstream brzegOdNas("PrzyNasPot", std::ios::binary);
	std::ofstream brzegDaleko("DalszyPot", std::ios::binary);

	for(int i=0; i<xElem; i++){
		for(int j=0; j<yElem; j++){
			brzegDol  << pot[i][j][0] << " ";
			brzegGora << pot[i][j][zEdge] << " ";
		}
		brzegDol  << std::endl;
		brzegGora << std::endl;
	}

	for(int i=0; i<xElem; i++){
		for(int k=0; k<zElem; k++){
			brzegPrawy 	<< pot[i][yEdge][k] 	<< " ";
			brzegLewy  << pot[i][0][k] 		<< " ";
			brzegDaleko << pot[xEdge][i][k] 	<< " ";
			brzegOdNas 	<< pot[0][i][k] 		<< " ";
		}
		brzegLewy 	<< std::endl;
		brzegPrawy 	<< std::endl;
		brzegOdNas  << std::endl;
		brzegDaleko << std::endl;
	}

	std::ofstream dens("koncentracjaPlanet", std::ios::binary);
	std::ofstream potFile("potencjal", std::ios::binary);

	for(int i=0; i<xElem; i++){
		for(int j=0; j<yElem; j++){
			for(int k=0; k<zElem; k++){
				dens		<< density[i][j][k] 	<< " ";
				potFile  << pot[i][j][k] 		<< " ";
			}
			dens 		<< std::endl;
			potFile 	<< std::endl;
		}
		dens << "\t";
		potFile << "\t";
	}

	potFile.close();
	dens.close();
	brzegDol.close();
	brzegGora.close();
	brzegDaleko.close();
	brzegOdNas.close();
	brzegLewy.close();
	brzegPrawy.close();
}

// TODO niedokładności na poziomie ~5% przy przypisaniu wartości początkowych, ~1/6 dx
