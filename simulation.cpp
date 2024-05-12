#include <iostream>
#include <fstream>
#include <cmath>

typedef long double ldouble;

template<typename T> T** alloc(int a, int b){
    T** temp = new T*[a];
    for(int i=0; i<a; i++){
        temp[i] = new T[b];
        for(int j=0; j<b; j++)
            temp[i][j] = 0;
    }
    return temp;
}

template<typename T> T*** alloc(int a, int b, int c){
    T*** temp = new T**[a];
    for(int i=0; i<a; i++){
        temp[i] = new T*[b];
        for(int j=0; j<b; j++){
            temp[i][j] = new T[c];
            for(int k=0; k<c; k++){
                temp[i][j][k] = 0;
            }
        }
    }
    return temp;
}

class wektor3d{
    public:
        long double x, y, z;
        wektor3d():x{0.0L}, y{0.0L}, z{0.0L}{}
        wektor3d(long double a, long double b, long double c):x{a}, y{b}, z{c}{}
        long double norm(){   return sqrt(x*x+y*y+z*z); }
        long double square2D(){ return x*x+y*y; }
        long double square(){ return x*x + y*y + z*z;}
        long double dot(wektor3d vec){ return x*vec.x + y*vec.y + z*vec.z; }
        wektor3d cross(wektor3d vec){ return wektor3d(y*vec.z-vec.y*z, vec.x*z-vec.z*x, x*vec.y-vec.x*y); }
        wektor3d operator*(wektor3d vec){ return wektor3d(vec.x*x, vec.y*y, vec.z*z);}
        wektor3d operator+(wektor3d vec){ return wektor3d(vec.x+x, vec.y+y, vec.z+z);}
        wektor3d operator-(wektor3d vec){ return wektor3d(x-vec.x, y-vec.y, z-vec.z);}
        wektor3d operator-(){ return wektor3d(-x, -y, -z); }
        template <typename T> wektor3d operator/(T a){
            return wektor3d(x/(long double)a, y/(long double)a, z/(long double)a);
        }
        template <typename T> wektor3d operator*(T a){
            return wektor3d(x*(long double)a, y*(long double)a, z*(long double)a);
        }
        void operator+=(wektor3d vec){ x+=vec.x; y+=vec.y; z+=vec.z; }
        template <typename T> friend wektor3d operator*(T a, wektor3d vec){
            return wektor3d(vec.x*(long double)a, vec.y*(long double)a, vec.z*(long double)a);
        }
        void operator=(wektor3d vec){ x=vec.x; y=vec.y; z=vec.z;}
};

inline void indices(int*, wektor3d);

// stałe
const ldouble G = 6.67e-11;
const ldouble M = 4.2e41;       // masa galaktyki
const ldouble D = 15*3.086e19;  // promień galaktyki

// parametry
const ldouble dt = 1e12;
const int iters = 3000;
const int N = 20000;
const int nx = 50;
const int nz = 1;

const int coIleResz = iters/100;
const int coIlePot  = iters/30;

// pomocnicze wielkości i stałe fizyczne
const ldouble dx = D/(nx+0.5);
const ldouble dwaDx = 2*dx;
const ldouble m = M/N;
const ldouble Gm = G*m;
const ldouble dt2 = dt*dt;
const ldouble Ckinet = m/2/dt/dt;
const int Nx = 2*nx+1;
const int Nz = 2*nz+1;
const int xEdge = Nx-1;
const int zEdge = Nz-1;

int main(){
    std::ofstream info("info",       std::ios::binary);
    std::ofstream loc( "Czatki",     std::ios::binary);
    std::ofstream pot( "Potencjal",  std::ios::binary);
    std::ofstream den( "gestosc",    std::ios::binary);
    std::ofstream ltv( "energie",    std::ios::binary);

    ldouble ek, ep, L, temp;

    info << dt << " " << Nx << " " << Nz << " " << dx << " " << iters << " " << coIleResz << " " << coIlePot << " " << N << std::endl; // zapisanie parametrów żeby odczytać je potem w pythonie

    int*** ro = alloc<int>(Nx, Nx, Nz);
    int** indx = alloc<int>(N, 3);
    bool* czyWsiatce = new bool[N];
    ldouble*** fi = alloc<ldouble>(Nx,Nx,Nz);
    ldouble*** tablica = alloc<ldouble>(Nx, Nx, Nz);
    wektor3d* rs = new wektor3d[N];
    wektor3d* rn = new wektor3d[N];
    wektor3d* rt = new wektor3d[N];
    wektor3d* a  = new wektor3d[N];

    for(int x=0; x<Nx; x++){ // tablicowanie odległości
        for(int y=0; y<Nx; y++){
            for(int z=0; z<Nz; z++){
                tablica[x][y][z] = Gm/(wektor3d(x, y, z)*dx).norm();
            }
        }
    }

    int x, y;
    for(int i=0; i<N/2; i++){ // pętla przypisująca położenia początkowe i indeksy, w których znajdują się cząstki
        do{
            x = (rand()%(xEdge-1) - nx+1), y = (rand()%(xEdge-1) - nx+1);
        }
        while( std::sqrt(std::pow(x, 2) + std::pow(y, 2)) > nx/2 ); 	// warunek dot. geometrii początkowej
        
        rs[2*i] = wektor3d(x+rand()/RAND_MAX-1/2, y+rand()/RAND_MAX-1/2, rand()%(zEdge-1) - nz + 1)*dx;
        rs[2*i+1] = -rs[2*i];

        indices(indx[2*i], rs[2*i]);
        indices(indx[2*i+1], rs[2*i+1]);
    }

    for(int i=0; i<N; i++){
        czyWsiatce[i] = true;
    }

    for(int i=0; i<N; i++){ 						// pętla odpowiedzialna za przesunięcie początkowe.
        ldouble rtemp = rs[i].norm()-0.75*dx; 				// to zmniejszać bo cząstki wypadają z siatki, zmniejszać czy zwiększać?
        int particles = 0;
        for(int j=0; j<N; j++){
            if(rs[j].norm()<rtemp && i!=j){
                particles += 1;
            }
        }
        ldouble w = 0.0L;
        if(rs[i].norm() > 1e15){
            w = sqrt(Gm*particles/pow(rs[i].norm(), 3));
        }
        else{
            w = 0.0L;   						// gdy r bliskie 0 to omega jest zerowana żeby jej do inf nie wybiło
        }
        rt[i] = wektor3d(rs[i].x*cos(dt*w) - rs[i].y*sin(w*dt), rs[i].y*cos(w*dt) + rs[i].x*sin(w*dt), rs[i].z); // obrót
    }

    for(int step=0; step<iters; step++){ 				// Pętla głównej symulacji
        for(int i=0; i<N; i++){ 					// znalezienie indeksów cząstek, przypisanie gęstości i sprawdzenie czy cząstki są w siatce.
            indices(indx[i], rt[i]);
            if(indx[i][0]>0 && indx[i][0]<xEdge && indx[i][1]>0 && indx[i][1]<xEdge && indx[i][2]<zEdge && indx[i][2]>0){
                czyWsiatce[i] = true;
                ro[indx[i][0]][indx[i][1]][indx[i][2]] += 1;
            }
            else{
                czyWsiatce[i] = false;
            }
        }

        for(int x=0; x<Nx; x++){ 					// potencjał na brzegu dla z = 0 i z = max
            for(int y=0; y<Nx; y++){
                fi[x][y][0] = 0.0L;
                fi[x][y][zEdge] = 0.0L;
                for(int i=1; i<x; i++){
                    for(int k=1; k<zEdge; k++){    
                        for(int j=1; j<y; j++){
                            temp = ro[i][j][k] * tablica[x-i][y-j][k];
                            fi[x][y][0] -= temp;
                            fi[x][y][zEdge] -= temp;
                        }
                        for(int j=y; j<xEdge; j++){
                            temp = ro[i][j][k] * tablica[x-i][j-y][k];
                            fi[x][y][0] -= temp;
                            fi[x][y][zEdge] -= temp;
                        }
                    }
                }
                for(int k=1; k<zEdge; k++){
                    for(int i=x; i<xEdge; i++){
                        for(int j=1; j<y; j++){
                            temp = ro[i][j][k] * tablica[i-x][y-j][k];
                            fi[x][y][0] -= temp;
                            fi[x][y][zEdge] -= temp;
                        }
                        for(int j=y; j<xEdge; j++){
                            temp = ro[i][j][k] * tablica[i-x][j-y][k];
                            fi[x][y][0] -= temp;
                            fi[x][y][zEdge] -= temp;
                        }
                    }
                }
                // for(int i=0; i<N; i++){ // dla małej ilość cząstek ta pętla jest szybsza.
                //     if(czyWsiatce[i]){
                //         temp = tablica[std::abs(x-indx[i][0])][std::abs(y-indx[i][1])][indx[i][2]];
                //         fi[x][y][0] -=  temp;
                //         fi[x][y][zEdge] -= temp;
                //     }
                // }
            }
        }

        for(int z=1; z<zEdge; z++){ 					// potencjał na brzegu dla y = 0 i y = max
            for(int x=0; x<Nx; x++){
                fi[x][0][z] = 0.0L;
                fi[x][xEdge][z] = 0.0L;
                for(int i=1; i<x; i++){
                    for(int j=1; j<xEdge; j++){
                        for(int k=1; k<z; k++){
                            temp = ro[i][j][k] * tablica[x-i][j][z-k];
                            fi[x][0][z] -= temp;
                            fi[x][xEdge][z] -= temp;
                        }
                        for(int k=z; k<zEdge; k++){
                            temp = ro[i][j][k] * tablica[x-i][j][k-z];
                            fi[x][0][z] -= temp;
                            fi[x][xEdge][z] -= temp;
                        }
                    }
                }
                for(int i=x; i<xEdge; i++){
                    for(int j=1; j<xEdge; j++){
                        for(int k=1; k<z; k++){
                            temp = ro[i][j][k] * tablica[i-x][j][z-k];
                            fi[x][0][z] -= temp;
                            fi[x][xEdge][z] -= temp;
                        }
                        for(int k=z; k<zEdge; k++){
                            temp = ro[i][j][k] * tablica[i-x][j][k-z];
                            fi[x][0][z] -= temp;
                            fi[x][xEdge][z] -= temp;
                        }
                    }
                }
                // for(int i=0; i<N; i++){
                //     if(czyWsiatce[i]){
                //         temp = tablica[std::abs(x-indx[i][0])][indx[i][1]][std::abs(z-indx[i][2])];
                //         fi[x][0][z] -= temp;
                //         fi[x][xEdge][z] -= temp;
                //     }
                // }
            }
        }

        for(int z=1; z<zEdge; z++){ 					// potencjał na brzegu na x = 0 i x = max
            for(int y=1; y<xEdge; y++){
                fi[0][y][z] = 0.0L;
                fi[xEdge][y][z] = 0.0L;
                for(int i=1; i<xEdge; i++){
                    for(int j=1; j<y; j++){
                        for(int k=1; k<z; k++){
                            temp = ro[i][j][k] * tablica[i][y-j][z-k];
                            fi[0][y][z] -= temp;
                            fi[xEdge][y][z] -= temp;
                        }
                        for(int k=z; k<zEdge; k++){
                            temp = ro[i][j][k] * tablica[i][y-j][k-z];
                            fi[0][y][z] -= temp;
                            fi[xEdge][y][z] -= temp;
                        }
                    }
                }
                for(int i=1; i<xEdge; i++){
                    for(int j=y; j<xEdge; j++){
                        for(int k=1; k<z; k++){
                            temp = ro[i][j][k] * tablica[i][j-y][z-k];
                            fi[0][y][z] -= temp;
                            fi[xEdge][y][z] -= temp;
                        }
                        for(int k=z; k<zEdge; k++){
                            temp = ro[i][j][k] * tablica[i][j-y][k-z];
                            fi[0][y][z] -= temp;
                            fi[xEdge][y][z] -= temp;
                        }
                    }
                }
                // for(int i=0; i<N; i++){
                //     if(czyWsiatce[i]){
                //         temp = tablica[indx[i][0]][std::abs(y-indx[i][1])][std::abs(z-indx[i][2])];
                //         fi[0][y][z] -= temp;
                //         fi[xEdge][y][z] -= temp;
                //     }
                // }
            }
        }

        for(int itd=0; itd<20; itd++){ 					// nadrelaksacja dla siatki poza brzegiem
            for(int x=1; x<xEdge; x++){
                for(int y=1; y<xEdge; y++){
                    for(int z=1; z<zEdge; z++){
                        fi[x][y][z] = -0.8*fi[x][y][z] + (1+0.8)/6.0*(fi[x+1][y][z]+fi[x-1][y][z]+fi[x][y+1][z]+fi[x][y-1][z]+fi[x][y][z+1]+fi[x][y][z-1]) - 1.8*4.0/6.0*M_PI*Gm/dx*ro[x][y][z];
                    }
                }
            }
        }

        if(step%coIlePot == 0){ 					// część zapisująca gęstość i potencjał
            for(int i=0; i<Nx; i++){
                for(int j=0; j<Nx; j++){
                    for(int k=0; k<Nz; k++){
                        pot << fi[i][j][k] << ' ';
                        den << ro[i][j][k] << ' ';
                    }
                    pot << '\n';
                    den << '\n';
                }
                den << '\t';
                pot << '\t';
            }
            den << '\r';
            pot << '\r';
        }

        for(int i=0; i<N; i++){ 					// tzw. mover
            if(czyWsiatce[i]){
                a[i].x = -(fi[indx[i][0]+1][indx[i][1]][indx[i][2]] - fi[indx[i][0]-1][indx[i][1]][indx[i][2]])/dwaDx;
                a[i].y = -(fi[indx[i][0]][indx[i][1]+1][indx[i][2]] - fi[indx[i][0]][indx[i][1]-1][indx[i][2]])/dwaDx;
                ro[indx[i][0]][indx[i][1]][indx[i][2]] = 0;
            }
            else{   							// siła dośrodkowa gdy jest poza siatką
                a[i] = -rt[i]*G*M/pow(rt[i].norm(), 3.0);
            }

            rn[i] = rt[i]*2 - rs[i] + a[i]*dt2;
            rs[i] = rt[i];
            rt[i] = rn[i];
        }

        if(step%coIleResz == 0){ 					// zapisanie położeń cząstek, energii i momentu pędu
            std::cout << (rt[0].dot(rt[0]-rs[0]))/(rt[0].norm()*(rt[0]-rs[0]).norm()) << std::endl;
            printf("Zrobiono: %.1lf %%, coIleResz = %d, step = %d\n", double(step*100)/iters, coIleResz, step);
            
            for(int i=0; i<N; i++){
                loc << rt[i].x << " " << rt[i].y << " " << a[i].x << " " << a[i].y << std::endl;
            }
            ep = 0.0, ek = 0.0, L = 0.0;

            for(int i=0; i<N; i++){
                ek += Ckinet * (rt[i]-rs[i]).square();
                if(czyWsiatce[i]){
                    ep += m * fi[indx[i][0]][indx[i][1]][indx[i][2]];
                }
                else{
                    ep -= Gm/rt[i].norm();
                }
                L  += m * (rt[i].cross(rt[i]-rs[i])).norm()/dt;
            }
            loc << "\t";
            ltv << ek << " " << ep << " " << L << std::endl;
        }
   }

   pot.close();
   den.close();
   loc.close();
   ltv.close();

   return 0;
}

inline void indices(int* indx, wektor3d r){ 				// funkcja licząca indeksy
		indx[0] = lround(r.x/dx) + nx;
		indx[1] = lround(r.y/dx) + nx;
		indx[2] = lround(r.z/dx) + nz;
}


// warunki początkowe w przypadku 2 cząstek
// rs[0].x = -6.0*dx;
// rs[1].x =  6.0*dx;

// for(int i=0; i<N; i++){
//     ldouble w = sqrt(G*m/(4*pow(rs[i].norm(), 3.0)));
//     rt[i] = wektor3d(rs[i].x*cos(dt*w) - rs[i].y*sin(w*dt), rs[i].y*cos(w*dt) + rs[i].x*sin(w*dt), 0); // o tyle obrót teoretycznie leci.
// }
