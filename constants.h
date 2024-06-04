#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>
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

#endif //CONSTANTS_H