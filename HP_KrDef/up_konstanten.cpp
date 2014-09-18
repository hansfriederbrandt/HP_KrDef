#include "stdafx.h"
#include "up_konstanten.h"

int P_Anz = 6;
double pi = M_PI;
double D = 75.15;       //[mm]; KR - Nenndurchmesser
double s = 0.210;        //[mm]; Spalthöhe
double ss = s / D;       //[]; Verhältnis von Spalthöhe zu Durchmesser
double L = pi*D - s;     //[mm]; Umfangslänge
double l = 0;          //[mm]; Elementlänge
double A = 0;        //[mm ^ 2] Querschnittsfläche
int n = 0;



//technische und physikalische Eigenschaften
double	 It = 8.374582659866258;           //[mm ^ 4] Flächentorsionsmoment
double	 E = 210000;        //[N / mm ^ 2]; E - Modul für Stahl
double	 G = 80000;         //[N / mm ^ 2]; E - Modul für Stahl
double	 rho = 7.85e-6;    //[kg / mm ^ 3]; Dichte für Stahl
double	 Jx = 0;     //[mm ^ 4] Querschnitts - Flächenträgheitsmoment in X - Richtung
double	 Jy = 0;     //[mm ^ 4] Querschnitts - Flächenträgheitsmoment in Y - Richtung
double	 EJ = 0;          //[Nmm ^ 2]
double	 EA = 0;           //[N]
double	 rhoA = 0;     //[kg / mm]
double   delta_XY = 0;

M_struct M1;


