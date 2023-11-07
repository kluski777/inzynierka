#ifndef _FUN_POMOCNICZE_
#define _FUN_POMOCNICZE_

typedef long double ldouble;
const int zElem = 4;
const int xElem = zElem*10, yElem = zElem*10;
   // założenie też jest takie że wymiary xElem == yElem
   // do tego liczba komórek jest zawsze parzysta
const int xEdge = xElem-1, yEdge = yElem-1, zEdge = zElem-1; 
const int N = xElem*zElem*yElem*10;

const int xElemD2 = xElem/2, yElemD2 = yElem/2, zElemD2 = zElem/2; // 20, 20, 2


///////////////////////////// Dane wielkości
const ldouble T = 6.93792e15;             // okres obrotu w sekundach
const ldouble M = 6e42;                   // Masa całej galaktyki
const ldouble mMin = M/N;                 // masa 1 cząstki w kg
const ldouble G = 6.67e-11;               // stała grawitacyjna w SI.
const ldouble w = 2*M_PI/T;               // częstość kołowa rad.
const ldouble dt = T/1000;                // krok czasowy, minimum 1000 kroków na obrót
const ldouble dx = 9.258e19/zElem;        // wymiar jednej kratki


//////////////////////////// Pomocnicze żeby przyspieszyć obliczenia
const ldouble dx2 = dx*dx;
const ldouble dx12 = .5*dx;
const ldouble coswt = cos(w*dt);
const ldouble sinwt = sin(w*dt);
const ldouble Gm = G*mMin;
const ldouble dt2dx = dt*dt/(2*dx); // rzędu 10^-32
const ldouble GM = G*M;
const ldouble m2 = mMin/2;
const ldouble dt2 = dt*dt;
const ldouble dt4 = 1.0/(4*dt2);
const ldouble potConst = (2.0/3.0)*M_PI*dx2*G*mMin/pow(dx,3.0);
const ldouble one6th = 1.0/6.0;


inline int xCoord2indx(ldouble coord){ return round((coord-dx12)/dx + xElemD2); }
inline int yCoord2indx(ldouble coord){ return round((coord-dx12)/dx + yElemD2); }
inline int zCoord2indx(ldouble coord){ return round((coord-dx12)/dx + zElemD2); }
inline ldouble xIndx2Coord(int indx){  return (indx-xElemD2)*dx + dx12; }
inline ldouble yIndx2Coord(int indx){  return (indx-yElemD2)*dx + dx12; }
inline ldouble zIndx2Coord(int indx){  return (indx-zElemD2)*dx + dx12; }

#endif

/*
std::ofstream potek("Potek.txt");
for(int wspX=0; wspX<xElem; wspX++){
for(int wspY=0; wspY<yElem; wspY++)for(int wspZ=0; wspZ<zElem; wspZ++)
	 potek << (20==wspX && 20==wspY ? 0 : Gm*density[wspX][wspY][1]*space[abs(20-wspX)][abs(20-wspY)][1]) << " ";
}
potek << std::endl;
}*/

/*                                                                                                           std::ofstream out1("Brzegowo.txt");
for(int i=0; i<xElem; i++){
for(int j=0; j<yElem; j++){
 out1 << pot[i][j][1] << " ";
}out1 << std::endl;}                                                                                                                                                                                                                         for(int i=0; i<xElem; i++){
for(int j=0; j<yElem; j++){
 ldouble r = sqrt(pow(xIndx2Coord(i), 2.0) + pow(yIndx2Coord(j), 2.0) + pow(zIndx2Coord(1), 2.));
 pot[i][j][1] = GM/(r);
}                                                                                                         }                                                                                                                                                                                                                         std::ofstream teo("GMr", std::ios::binary);
for(int i=0; i<xElem; i++){
for(int j=0; j<yElem; j++)
 teo << pot[i][j][1] << " ";
teo << std::endl;                                                                                         }*/

   /*
   std::ofstream out("filik.txt");
   for(int i=1; i<xElem-1; i++){
      ldouble L = (pot[i-1][10][1]+pot[i+1][10][1]+pot[i][11][1]+pot[i][9][1]+pot[i][10][0]+pot[i][10][2]-6*pot[i][10][1])/dx2;
      ldouble P = 4*M_PI*G*mMin/pow(dx,3.0)*density[i][10][1];
      out << i << " " << P << " " << L << std::endl;
   }*/
