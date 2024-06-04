#include <iostream>
#include "animation.hpp"
#include "constants.hpp"  

using std::cout;
using std::endl;

void StartAnimation(void) {
    cout << "set terminal gif animate" << endl;
    cout << "set output 'phases.gif'" << endl;
    cout << "unset key" << endl;
    cout << "set size ratio 4" << endl;
    cout << "set xrange [0:" << Lx << "]" << endl;
    cout << "set yrange [0:" << Ly << "]" << endl;
    cout << "set cbrange["<<rho_l-0.1*rho_l<<":"<<rho_h+0.1*rho_h<<"]" << endl;
}

void StartFrame(void) {
    cout << "plot 0,0 ";
    cout<<", 'data.dat' w image";
}

void EndFrame(void) {
    cout << endl;
}