#ifndef _FUN_POMOCNICZE_
#define _FUN_POMOCNICZE_

/* Sprawdzić brzegi i dlaczego cząstki przyspieszają?
 * 
 * 
 * 
 * 
*/

typedef long double ldouble; 	// zmienić zmienną na jakąś dokładniejszą.
const int zElem = 8; 			// zawsze parzysta ma być!
const int xElem = zElem*10, yElem = zElem*10;

const int xEdge = xElem-1, yEdge = yElem-1, zEdge = zElem-1; 
const int N = xElem*zElem*yElem*10;

const int xElemD2 = xElem/2, yElemD2 = yElem/2, zElemD2 = zElem/2;
const int xElemD4 = xElemD2/2, yElemD4 = yElemD2/2, zElemD4 = zElemD2/2;

///////////////////////////// Dane wielkości
const ldouble okres = 6.93792e15;			// okres obrotu w sekundach
const ldouble M = 2.9835e42;              // Masa całej galaktyki
const ldouble mMin = M/N;                 // masa 1 cząstki w kg
const ldouble G = 6.6743e-11; 				// stała grawitacyjna w SI.
const ldouble dt = okres/1000;           // krok czasowy, minimum 1000 kroków na obrót
const ldouble dx = 9.25703274384e19/zElem;// wymiar jednej kratki


//////////////////////////// Pomocnicze żeby przyspieszyć obliczenia
const ldouble dx2 = dx*dx;
const ldouble dx12 = .5*dx;
//const ldouble coswt = cos(w*dt);
//const ldouble sinwt = sin(w*dt);
const ldouble Gm = G*mMin;
const ldouble dt2dx = dt*dt/(2*dx); // rzędu 10^-32
const ldouble GM = G*M;
const ldouble m2 = mMin/2;
const ldouble dt2 = dt*dt;
const ldouble dt4 = 1.0/(4.0*dt2);
const ldouble potConst = (2.0/3.0)*M_PI*dx2*G*mMin/pow(dx,3.0);
const ldouble one6th = 1.0/6.0;


inline int xCoord2indx(ldouble coord){ return round((coord-dx12)/dx + xElemD2); }
inline int yCoord2indx(ldouble coord){ return round((coord-dx12)/dx + yElemD2); }
inline int zCoord2indx(ldouble coord){ return round((coord-dx12)/dx + zElemD2); }
inline ldouble xIndx2Coord(int indx){  return (indx-xElemD2+.5)*dx; }
inline ldouble yIndx2Coord(int indx){  return (indx-yElemD2+.5)*dx; }
inline ldouble zIndx2Coord(int indx){  return (indx-zElemD2+.5)*dx; }
inline float indxDstce2(int i, int j, int k){	return pow(i-xElemD2+.5, 2) + pow(j-yElemD2+.5, 2) + pow(k-zElemD2+.5, 2); }

template <typename T> void saveMap(T*** arr, int k, std::string title){
	std::ofstream outF(title, std::ios::binary);
	T sum=0;
	for(int i=0; i<xElem; i++){
		for(int j=0; j<yElem; j++){
			outF << arr[i][j][k] << " ";
			sum += arr[i][j][k];
		}
		outF << std::endl;
	}
	printf("Suma: %d\n", sum);
	outF.close();
}

#endif
