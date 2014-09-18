#ifndef UP_VARIABLEN_H
#define UP_VARIABLEN_H

#include <math.h>
#include <Eigen/Eigen>
using namespace Eigen;

//Fl�chenberechnung/Bearbeitung des Querschnitts
extern int P_Anz;

//geometrische Gr��en
extern double pi;
extern double D;       //[mm]; KR - Nenndurchmesser
extern double s;        //[mm]; Spalth�he
extern double ss;       //[]; Verh�ltnis von Spalth�he zu Durchmesser
extern double L;     //[mm]; Umfangsl�nge
extern double l_Balk; //[mm] L�nge der einzelnen Balkenelemente f�r bestes Verh�ltnis f�r FEM    %n = input('Balkenanzahl n = ');   %Anzahl der Balkenelemente in die Ring eingeteilt wird.
extern int n; //Anzahl der Balkenelemente
extern double l;          //[mm]; Elementl�nge
extern double A;        //[mm ^ 2] Querschnittsfl�che

//technische und physikalische Eigenschaften
extern double	 It;           //[mm ^ 4] Fl�chentorsionsmoment
extern double	 Jx;     //[mm ^ 4] Querschnitts - Fl�chentr�gheitsmoment in X - Richtung
extern double	 Jy;     //[mm ^ 4] Querschnitts - Fl�chentr�gheitsmoment in Y - Richtung
extern double	 E;        //[N / mm ^ 2]; E - Modul f�r Stahl
extern double	 G;         //[N / mm ^ 2]; E - Modul f�r Stahl
extern double	 EJ;          //[Nmm ^ 2]
extern double	 EA;           //[N]
extern double	 rho;    //[kg / mm ^ 3]; Dichte f�r Stahl
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

