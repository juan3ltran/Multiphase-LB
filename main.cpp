
#include"animation.h"
#include"Latticeboltzmann.h"


int main(void){
    StartAnimation();
    LatticeBoltzmann instability;
    int t, tmax = 400000;
    double Ux0 = 0; double Uy0 = 0;
    instability.Start(Ux0, Uy0);
    for(t=0;t<tmax;t++){
       instability.Collision();
       instability.ImposeFields();
       instability.Advection();
       if (t%500 == 0){
    
       instability.print("data.dat");
       StartFrame();
       EndFrame();}
    }
    
    return 0;
}