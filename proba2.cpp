#include <iostream>
#include <cmath>
#include <fstream>

int main() {
    const int nx = 12;
    const int nz = 1;
    const int n = 2;
    const double rsize = 30 * 3.086e19;
    const double g = 6.67e-11;
    const double dx = rsize / Nx;
    const double xmc = 6e42 / n;
    const double pi = 4 * atan(1.0);
    const double dt = 1e10;

	 std::ofstream outputFile("Czastki");
    
	 double phi[2 * nx + 1][2 * nx + 1][2 * nz + 1];
    double ro[2 * nx + 1][2 * nx + 1][2 * nz + 1];
    double rs[n][3], a[n][3], v[n][3], rn[n][3], r[n][3], nk[n][3];

    // Initialize rs, v, and r
    rs[0][0] = -3.5 * dx;
    rs[0][1] = 0;
    rs[1][0] = 3.5 * dx;
    rs[1][1] = 0;

    v[0][1] = sqrt(g * xmc / 7 / dx) / 1.41;
    v[1][1] = -sqrt(g * xmc / 7 / dx) / 1.41;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < 3; ++j) {
            r[i][j] = rs[i][j] + v[i][j] * dt;
        }
    }

    // Initialize phi
    for (int i = -nx; i <= nx; ++i) {
        for (int j = -nx; j <= nx; ++j) {
            for (int k = -nz; k <= nz; ++k) {
                phi[i + nx][j + nx][k + nz] = 0;
            }
        }
    }


    for (int iter = 1; iter <= 100000; ++iter) {
        // Reset ro
        for (int i = 0; i < n; ++i) {
            for (int ix = -nx; ix <= nx; ++ix) {
                for (int iy = -nx; iy <= nx; ++iy) {
                    for (int iz = -nz; iz <= nz; ++iz) {
                        double x = ix * dx;
                        double y = iy * dx;
                        double z = iz * dx;
                        double xd = std::abs(r[i][0] - x);
                        double yd = std::abs(r[i][1] - y);
                        double zd = std::abs(r[i][2] - z);

                        if (xd < dx / 2.0 && yd < dx / 2.0 && zd < dx / 2.0) {
                            ro[ix + nx][iy + nx][iz + nz] += xmc / std::pow(dx, 3);
                        }
                    }
                }
            }
        }

		for (int ix = -nx; ix <= nx; ++ix) {
    		for (int iy = -nx; iy <= nx; ++iy) {
        		for (int iz = -nz; iz <= nz; ++iz) {
           		 if (std::abs(iz) == nz || std::abs(ix) == nx || std::abs(iy) == nx) {
               	 phi[ix + nx][iy + nx][iz + nz] = 0;
                		for (int i = 0; i < n; ++i) {
                    		double rc = std::pow(ix * dx - r[i][0], 2) + std::pow(iy * dx - r[i][1], 2) + std::pow(iz * dx - r[i][2], 2);
                    		rc = std::sqrt(rc);
                    		phi[ix + nx][iy + nx][iz + nz] -= g * xmc / rc;
               		 }
                		// Additional write statement (commented out in Fortran)
                	// std::cout << ix << " " << iy << " " << iz << " " << phi[ix + nx][iy + nx][iz + nz] << std::endl;
           		}
        		}
    		}
		}

	for (int itp = 1; itp <= 200; ++itp) {
   	 for (int ix = -nx + 1; ix <= nx - 1; ++ix) {
      	  for (int iy = -nx + 1; iy <= nx - 1; ++iy) {
         	   	for (int iz = -nz + 1; iz <= nz - 1; ++iz) {
              	  phi[ix + nx][iy + nx][iz + nz] = phi[ix + nx][iy + nx][iz + nz] * (-0.8) +
                  	  (1 + 0.8) * (1.0 / 6 * (phi[ix - 1 + nx][iy + nx][iz + nz] +
                           	                 phi[ix + 1 + nx][iy + nx][iz + nz] +
                        	                    phi[ix + nx][iy + 1 + nx][iz + nz] +
                     	                       phi[ix + nx][iy - 1 + nx][iz + nz] +
                  	                          phi[ix + nx][iy + nx][iz + 1 + nz] +
               	                             phi[ix + nx][iy + nx][iz - 1 + nz]) -
            	                     4 * pi / 6 * std::pow(dx, 2) * g * ro[ix + nx][iy + nx][iz + nz]);
         	   }
      	  }
   	 }
	}

for (int i = 0; i < n; ++i) {
    for (int ix = -nx; ix <= nx; ++ix) {
        for (int iy = -nx; iy <= nx; ++iy) {
            for (int iz = -nz; iz <= nz; ++iz) {
                double x = ix * dx;
                double y = iy * dx;
                double z = iz * dx;
                double xd = std::abs(r[i][0] - x);
                double yd = std::abs(r[i][1] - y);
                double zd = std::abs(r[i][2] - z);

                if (xd < dx / 2.0 && yd < dx / 2.0 && zd < dx / 2.0) {
                    nk[i][0] = ix + nx;
                    nk[i][1] = iy + nx;
                    nk[i][2] = iz + nz;
                    // Additional write statement (commented out in Fortran)
                    // std::cout << i + 1 << " " << ix << " " << iy << " " << iz << std::endl;
                }
            }
        }
    }
}

for (int i = 0; i < n; ++i) {
    int ix = nk[i][0];
    int iy = nk[i][1];
    int iz = nk[i][2];
    a[i][0] = -(phi[ix + 1][iy][iz] - phi[ix - 1][iy][iz]) / (2 * dx);
    a[i][1] = -(phi[ix][iy + 1][iz] - phi[ix][iy - 1][iz]) / (2 * dx);
    a[i][2] = 0;
}

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < 3; ++j) {
                rn[i][j] = 2 * r[i][j] - rs[i][j] + a[i][j] * std::pow(dt, 2);
            }
        }
			
			for(int i=0; i < n; i++){
				for(int j=0; j < 3; j++){
        			rs[i][j] = r[i][j];
        			r[i][j] = rn[i][j];
				}
			}

        double ek = 0;
        double ep = 0;

        // Additional calculations...

        if (iter % 1000 == 0) {
            outputFile << r[0][0] << " " << r[0][1] << " " << r[1][0] << " " << r[1][1] << " " << ek << " " << ep << " " << ek + ep << std::endl;
        }
    }

    outputFile.close();
    
	 return 0;
}

