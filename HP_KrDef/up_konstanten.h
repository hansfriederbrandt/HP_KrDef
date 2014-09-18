#ifndef UP_VARIABLEN_H
#define UP_VARIABLEN_H

#include <math.h>
#include <Eigen/Eigen>
using namespace Eigen;

//Flächenberechnung/Bearbeitung des Querschnitts
extern int P_Anz;

//geometrische Größen
extern double pi;
extern double D;       //[mm]; KR - Nenndurchmesser
extern double s;        //[mm]; Spalthöhe
extern double ss;       //[]; Verhältnis von Spalthöhe zu Durchmesser
extern double L;     //[mm]; Umfangslänge
extern double l_Balk; //[mm] Länge der einzelnen Balkenelemente für bestes Verhältnis für FEM    %n = input('Balkenanzahl n = ');   %Anzahl der Balkenelemente in die Ring eingeteilt wird.
extern int n; //Anzahl der Balkenelemente
extern double l;          //[mm]; Elementlänge
extern double A;        //[mm ^ 2] Querschnittsfläche

//technische und physikalische Eigenschaften
extern double	 It;           //[mm ^ 4] Flächentorsionsmoment
extern double	 Jx;     //[mm ^ 4] Querschnitts - Flächenträgheitsmoment in X - Richtung
extern double	 Jy;     //[mm ^ 4] Querschnitts - Flächenträgheitsmoment in Y - Richtung
extern double	 E;        //[N / mm ^ 2]; E - Modul für Stahl
extern double	 G;         //[N / mm ^ 2]; E - Modul für Stahl
extern double	 EJ;          //[Nmm ^ 2]
extern double	 EA;           //[N]
extern double	 rho;    //[kg / mm ^ 3]; Dichte für Stahl
extern double	 rhoA;     //[kg / mm]

extern double delta_XY;

struct M_struct{
	MatrixXd KM_sys;
	MatrixXd MM_sys;
	MatrixXd GM_sys;
	MatrixXd TM_sys;
};

struct F{ double F_Res_x; };

#endif

