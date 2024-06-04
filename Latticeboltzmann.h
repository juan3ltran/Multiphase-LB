#include"constants.h"
#include"vector.h"

class LatticeBoltzmann{
private:
    double w[Q]; //weights
    int Vx[Q], Vy[Q]; //velocity vectors
    vector3D V[Q];
    
    //Distribution Functions
    double *f_bar, *f_bar_new;
    double *g_bar, *g_bar_new;

    //Potential
    double *PSI;


public:
    LatticeBoltzmann(void);
    ~LatticeBoltzmann(void);
    int n(int ix, int iy, int i){return (ix*Ly+iy)*Q+i;};
    double Psi(int ix, int iy);
    double rho(int ix, int iy, double phi);
    double phi(int ix, int iy, bool Usenew);
    double pressure(int ix, int iy, vector3D U, vector3D grad, bool UseNew);
    double feq(double rho0, vector3D U0, int i);
    double geq(double p0, double rho0, vector3D U0, int i);
    double gradient_x(int ix, int iy);
    double gradient_y(int ix, int iy);
    double Gamma(vector3D U0, int i);
    vector3D U(int ix, int iy, vector3D G, double rho0, bool UseNew);
    vector3D gradient(int ix, int iy);
    void Collision(void);
    void Advection(void);
    void ImposeFields(void);
    void Start(double Ux0, double Uy0);
    void print( const char * NombreArchivo);
    
};