// HP_KrDef.cpp : Definiert den Einstiegspunkt für die Konsolenanwendung.
//


#include "stdafx.h"
#include "up_konstanten.h"
#include "up_startberechnungen.h"
#include "up_matritzen.h"
#include "up_kraftvektor.h"
#include "up_ode.h"
#include "up_rk4_abm.h"




int main()
{
	double h = 1e-7;
	VectorXd F_Eingang(3);
	double F_grad;
	int j;
	F_Eingang(0) = 0;
	F_Eingang(1) = -113;
	F_Eingang(2) = (F_Eingang(1) - F_Eingang(0)) / h;



	M_struct M1;

	Vector2d PS;

	j = 1;


	up_startberechnungen(PS);
	up_matritzen(M1);
	MatrixXd Y(4 * (n + 1), 2*4);
	Y.fill(0);
	up_rk4_abm(F_Eingang, h, Y, M1, j);






}