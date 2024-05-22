#include <iostream>
#include <fstream>
#include <cmath>
#include "vector.h"
#include <omp.h>
using std::cout, std::endl;

const int Lx = 128;
const int Ly = 256;

const int Q = 9;


const double a = 6.93;//7.225;
const double RT = a/12.0;//.5;
const double c = std::sqrt(3.0*RT);
const double phi_l = 0.0317;//0.0;
const double phi_h = 0.2554;//0.278575;
const double rho_l = 1;
const double rho_h =  3;
const double kappa = 0;
const double g = 0.008/Ly;
const double ps = 0.018;



const double tau =  0.56;
const double Utau = 1.0/tau;
const double DtaumUsDtau = (2*tau - 1 )/(2*tau);

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

LatticeBoltzmann::LatticeBoltzmann(void){
    //set wheights
    w[0] = 4.0/9;
    w[1] = w[2] = w[3] = w[4] = 1.0/9;
    w[5] = w[6] = w[7] = w[8] =1.0/36;
    //set velocity vectors
    Vx[0] = 0; Vx[1] = 1; Vx[2] = 0; Vx[3] = -1; Vx[4] = 0;
    Vy[0] = 0; Vy[1] = 0; Vy[2] = 1; Vy[3] = 0; Vy[4] = -1;

    Vx[5] = 1; Vx[6] = -1; Vx[7] = -1; Vx[8] = 1;
    Vy[5] = 1; Vy[6] = 1; Vy[7] = -1; Vy[8] = -1; 
    
    for (int i = 0; i <Q; i++){
        V[i].load(Vx[i],Vy[i],0);
    }

    //Dynamic arrays
    int ArraySize = Lx*Ly*Q;
    f_bar = new double[ArraySize]; f_bar_new = new double[ArraySize];
    g_bar = new double[ArraySize]; g_bar_new = new double[ArraySize];
    PSI = new double [Lx*Ly];
}

LatticeBoltzmann::~LatticeBoltzmann(void){
    delete[] f_bar; delete[] f_bar_new;
    delete[] g_bar; delete[] g_bar_new;
    delete[] PSI;
}

double LatticeBoltzmann::Psi(int ix,int iy){
    double  phi0=LatticeBoltzmann::phi(ix, iy,false), phi2 = phi0*phi0;
    return phi2*RT*(4.0-2*phi0)/pow(1-phi0,3)-a*phi2;
}

double LatticeBoltzmann::rho(int ix, int iy, double phi){
    return rho_l + (phi-phi_l)*(rho_h-rho_l)/(phi_h-phi_l);
}

double LatticeBoltzmann::phi(int ix, int iy, bool UseNew){
    double sum ;int i,n0;
    for(sum=0, i =0; i<Q; i++){
        n0 = LatticeBoltzmann::n(ix, iy, i);
        if(UseNew){
            sum += f_bar_new[n0];
        }
        else{
            sum += f_bar[n0];
        }
    }
    return sum;
}

vector3D LatticeBoltzmann::U(int ix, int iy, vector3D G, double rho0, bool UseNew){
    vector3D sum ;int i,n0;
    for(sum.load(0,0,0), i =0; i<Q; i++){
        n0 = LatticeBoltzmann::n(ix, iy, i);
        if(UseNew){
            sum += V[i]*g_bar_new[n0];
        }
        else{
            sum += V[i]*g_bar[n0];
        }
    }
    return (sum + 0.5*RT*G)* 1/(rho0*RT);
}

double LatticeBoltzmann::pressure(int ix, int iy, vector3D U, vector3D grad, bool UseNew){
    double sum = 0; int i, n0;
    for(i =0; i<Q; i++){
        n0 = LatticeBoltzmann::n(ix, iy, i);
        if(UseNew){
            sum += g_bar_new[n0];
        }
        else{
            sum += g_bar[n0];
        }
    }
    return sum - 0.5*U*grad;
}


double LatticeBoltzmann::feq(double phi0, vector3D U0, int i){
    double UdotVi = U0*V[i], U2 = U0.norm2();
    return phi0*w[i]*(1.0 + (3.0*UdotVi)/(c*c) + (4.5*UdotVi*UdotVi)/(pow(c,4.0)) - 1.5*U2/(c*c));

}

double LatticeBoltzmann::geq(double p0, double rho0, vector3D U0, int i){
    double UdotVi = U0*V[i], U2 = U0.norm2();
    return w[i]*(p0 + rho0*((3.0*UdotVi)/(c*c) + (4.5*UdotVi*UdotVi)/(pow(c,4.0)) - 1.5*U2/(c*c)));
}

double LatticeBoltzmann::Gamma(vector3D U0, int i){
    double UdotVi = U0*V[i], U2 = U0.norm2();
    return w[i]*(1.0 + (3.0*UdotVi)/(c*c) + (4.5*UdotVi*UdotVi)/(pow(c,4.0)) - 1.5*U2/(c*c));
}


double LatticeBoltzmann::gradient_y(int ix, int iy){
    double sum_y = 0;
    for (int i= 0;i<9;i++){
        sum_y+= Vy[i]*(8.0*PSI[((ix+Vx[i] +Lx)%Lx)*Ly +((iy+Vy[i] + Ly)%Ly)] - PSI[((ix+2*Vx[i] +Lx)%Lx)*Ly +((iy+2*Vy[i] + Ly)%Ly)]);
    }
    int j1 = (iy + 1) % Ly; 
    int j2 = (iy + 2) % Ly;
    int j3 = (iy + 3) % Ly;
    
    int jm1 = (iy - 1 + Ly) % Ly; 
    int jm2 = (iy - 2 + Ly) % Ly;
    int jm3 = (iy - 3 + Ly) % Ly;

    double dy, grad_y;

    dy = 1.0/60.0 *( -PSI[ix*Ly+jm3] + 9.0*PSI[ix*Ly+jm2] - 45.0*PSI[ix*Ly+jm1] + 45.0*PSI[ix*Ly+j1] - 9.0*PSI[ix*Ly+j2] + PSI[ix*Ly+j3]);
    
    grad_y = 0.75*dy +  0.25*(1.0/(36.0*c))*sum_y;
    return grad_y;
}

double LatticeBoltzmann::gradient_x(int ix, int iy){
    double sum_x = 0;
    for (int i= 0;i<9;i++){
        sum_x += Vx[i]*(8.0*PSI[((ix+Vx[i] +Lx)%Lx)*Ly +((iy+Vy[i] + Ly)%Ly)] - PSI[((ix+2*Vx[i] +Lx)%Lx)*Ly +((iy+2*Vy[i] + Ly)%Ly)]);
    }
    

    int i1 = (ix + 1) % Lx; // siguiente fila, con envoltura
    int i2 = (ix + 2) % Lx;
    int i3 = (ix + 3) % Lx;

    int im1 = (ix - 1 + Lx) % Lx; // fila anterior, con envoltura
    int im2 = (ix - 2 + Lx) % Lx;
    int im3 =(ix - 3 + Lx) % Lx;

    double dx, grad_x;

    dx = 1.0/60.0 *( -PSI[im3*Ly+iy] + 9.0*PSI[im2*Ly+iy] - 45.0*PSI[im1*Ly+iy] + 45.0*PSI[i1*Ly+iy] - 9.0*PSI[i2*Ly+iy] + PSI[i3*Ly+iy]);
    grad_x = 0.75*dx + 0.25*(1.0/(36.0*c))*sum_x;
    return grad_x;
}

vector3D LatticeBoltzmann::gradient(int ix, int iy){
    double grad_x, grad_y;
    vector3D grad;
    grad_x = LatticeBoltzmann::gradient_x(ix,iy); grad_y = LatticeBoltzmann::gradient_y(ix,iy);
    grad.load(grad_x, grad_y,0);
    return grad;
}


void LatticeBoltzmann::Start(double Ux0, double Uy0){
    int ix, iy, i ,n0;
    double rho0, phi0, P0;
    vector3D U0;

    for(ix=0;ix<Lx;ix++){
        for(iy=0;iy<Ly;iy++){
            phi0 = (phi_h - phi_l)/M_PI * atan(double(iy)-Lx/20.0*cos(((2*M_PI)/(Lx))*double(ix))-double(Ly)/1.2) + (phi_h + phi_l)/2.0;
            rho0 = LatticeBoltzmann::rho(ix, iy, phi0);
            P0 = ps;
            U0.load(Ux0, Uy0,0);
            for(i=0;i<Q;i++){
                n0 = LatticeBoltzmann::n(ix,iy,i);
                f_bar[n0] = f_bar_new[n0] = LatticeBoltzmann::feq(phi0, U0, i);
                g_bar[n0] = g_bar_new[n0] =  LatticeBoltzmann::geq(P0, rho0, U0, i);
            }
        }
        
    }
}

void LatticeBoltzmann::Collision(void){
    
    int ix, iy, i, n0; 
    #pragma omp parallel for collapse(2)  
    for(ix=0;ix<Lx;ix++){
        for(iy=0;iy<Ly;iy++){   
            PSI[ix*Ly+iy] = LatticeBoltzmann::Psi(ix,iy);         
        }
    }

    double gammau, gamma0, phi0, rho0, P0;
    vector3D grad, G, U;
    #pragma omp parallel for collapse(2) private(phi0, grad, rho0, G, U, P0, i, n0, gamma0, gammau)
    for(int ix = 0; ix < Lx; ix++) {
        for(int iy = 0; iy < Ly; iy++) {  
            phi0 = LatticeBoltzmann::phi(ix, iy, false);
            grad = LatticeBoltzmann::gradient(ix, iy);
            rho0 = LatticeBoltzmann::rho(ix, iy, phi0);
            G.load(0, -rho0 * g, 0);
            U = LatticeBoltzmann::U(ix, iy, G, rho0, false); 
            P0 = LatticeBoltzmann::pressure(ix, iy, U, grad, false);

            for(int i = 0; i < Q; i++) {
                n0 = LatticeBoltzmann::n(ix, iy, i);
                double feq0 = LatticeBoltzmann::feq(phi0, U, i);
                double geq0 = LatticeBoltzmann::geq(P0, rho0, U, i);
                gamma0 = w[i];
                gammau = LatticeBoltzmann::Gamma(U, i);

                f_bar_new[n0] = f_bar[n0] - Utau * (f_bar[n0] - feq0) - DtaumUsDtau * ((V[i] - U) * (grad * gammau * (1 / RT)));
                g_bar_new[n0] = g_bar[n0] - Utau * (g_bar[n0] - geq0) + DtaumUsDtau * ((V[i] - U) * (gammau * G - (gammau - gamma0) * grad));
            }
        }
    }
}

void LatticeBoltzmann::ImposeFields(void){
    int i, ix, iy, n0; vector3D G, U;
    U.load(0,0,0);
    for (ix = 0; ix < Lx; ix++){
        for (iy = 0; iy < 1; iy ++){
            for(i=0;i<Q;i++){
                n0 = LatticeBoltzmann::n(ix,iy,i);
                f_bar_new[n0] = LatticeBoltzmann::feq(phi_l,U,i);
                g_bar_new[n0] = LatticeBoltzmann::geq(ps,rho_l,U,i);
            }
        }

        for (iy = Ly -1; iy < Ly; iy ++){
            for(i=0;i<Q;i++){
                n0 = LatticeBoltzmann::n(ix,iy,i);
                f_bar_new[n0] = LatticeBoltzmann::feq(phi_h,U,i);
                g_bar_new[n0] = LatticeBoltzmann::geq(ps,rho_h,U,i);
            }
        }
    }
    

}


void LatticeBoltzmann::Advection(void){
    int ix, iy, i ,n0, ixnext, iynext ,n0next;
    #pragma omp parallel for collapse(3) private(ixnext, iynext, n0)
    for(ix=0;ix<Lx;ix++){
        for(iy=0;iy<Ly;iy++){
            for(i=0;i<Q;i++){
                ixnext = (ix + Vx[i] + Lx)%Lx; iynext = (iy+Vy[i]+Ly)%Ly;
                n0 = LatticeBoltzmann::n(ix, iy, i); n0next = LatticeBoltzmann::n(ixnext, iynext, i);
                f_bar[n0next] = f_bar_new[n0];
                g_bar[n0next] = g_bar_new[n0];
            }
        }
    }
}

void LatticeBoltzmann::print(const char * Namefile){
    std::ofstream MyFile(Namefile); double rho0, phi, pressure; int ix, iy;
    for(ix=0;ix<Lx;ix+=1){
        for(iy=0;iy<Ly;iy+=1){
            vector3D U,G, grad;
            phi = LatticeBoltzmann::phi(ix,iy,true);
            rho0 = LatticeBoltzmann::rho(ix,iy,phi);
            
            MyFile<<ix<<" "<<iy<<" "<<rho0<<std::endl;
               
        }
        
        MyFile<<std::endl;
    }
    MyFile.close();
}

void StartAnimation(void){
  cout<<"set terminal gif animate"<<endl; 
  cout<<"set output 'fases.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set size ratio 4"<<endl;
  cout<<"set xrange [0:"<<Lx<<"]"<<endl;
  cout<<"set yrange [0:"<<Ly<<"]"<<endl;
  cout<<"set cbrange[0.9:3.1]"<<endl;
}
void StartFrame(void){
    cout<<"plot 0,0 ";
}
void EndFrame(void){
    cout<<endl;
}

int main(void){
    StartAnimation();
    LatticeBoltzmann instability;
    int t, tmax = 400000;
    double Ux0 = 0; double Uy0 = 0;
    instability.Start(Ux0, Uy0);
    for(t=0;t<tmax;t++){
        //std::cout<<std::endl<<"colision "<<t+1<<std::endl;
       instability.Collision();
       instability.ImposeFields();
       instability.Advection();
       if (t%500 == 0){
    
       instability.print("data.dat");
       StartFrame();
       cout<<", 'data.dat' w image";
       EndFrame();}
    }
    //inestability.print("datos.dat");
    return 0;
}
