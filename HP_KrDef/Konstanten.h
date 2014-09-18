//Konstanten.h
#define _USE_MATH_DEFINES
#include <math.h>
#include <Eigen/Eigen>
using namespace Eigen;

#if !defined(Konstanten)
#define Konstanten 1

static double pi = M_PI;
static double D = 75.15;       //[mm]; KR - Nenndurchmesser
static double s = 0.210;        //[mm]; Spalthöhe
double ss = s / D;       //[]; Verhältnis von Spalthöhe zu Durchmesser
static double L = pi*D - s;     //[mm]; Umfangslänge



//	 //technische und physikalische Eigenschaften
static double	 It = 8.374582659866258;           //[mm ^ 4] Flächentorsionsmoment
static double	 E = 210000;        //[N / mm ^ 2]; E - Modul für Stahl
static double	 G = 80000;         //[N / mm ^ 2]; E - Modul für Stahl
static double	 rho = 7.85e-6;    //[kg / mm ^ 3]; Dichte für Stahl



#endif