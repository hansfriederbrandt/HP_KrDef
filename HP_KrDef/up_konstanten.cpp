#include "stdafx.h"
#include "up_konstanten.h"

int P_Anz = 6;
double pi = M_PI;
double D = 75.15;       //[mm]; KR - Nenndurchmesser
double s = 0.210;        //[mm]; Spalth�he
double ss = s / D;       //[]; Verh�ltnis von Spalth�he zu Durchmesser
double L = pi*D - s;     //[mm]; Umfangsl�nge
double l = 0;          //[mm]; Elementl�nge
double A = 0;        //[mm ^ 2] Querschnittsfl�che
int n = 0;



//technische und physikalische Eigenschaften
double	 It = 8.374582659866258;           //[mm ^ 4] Fl�chentorsionsmoment
double	 E = 210000;        //[N / mm ^ 2]; E - Modul f�r Stahl
double	 G = 80000;         //[N / mm ^ 2]; E - Modul f�r Stahl
double	 rho = 7.85e-6;    //[kg / mm ^ 3]; Dichte f�r Stahl
double	 Jx = 0;     //[mm ^ 4] Querschnitts - Fl�chentr�gheitsmoment in X - Richtung
double	 Jy = 0;     //[mm ^ 4] Querschnitts - Fl�chentr�gheitsmoment in Y - Richtung
double	 EJ = 0;          //[Nmm ^ 2]
double	 EA = 0;           //[N]
double	 rhoA = 0;     //[kg / mm]
double   delta_XY = 0;

M_struct M1;


