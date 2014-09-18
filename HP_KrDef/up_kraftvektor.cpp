#include "stdafx.h"
#include "up_konstanten.h"
#include "up_kraftvektor.h"



VectorXd up_kraftvektor(double &F_Res, M_struct &M1){
	//function[F_sys] = UP_Kraftvektor(F_Res_x, T, n, l, delta_XY)
	//Kraftvektortransformation

	//Eingangsdaten
	//Fn-- > Gesamtkraft Fn auf Knoten verteilen.

	MatrixXd TMG;
	TMG = M1.TM_sys;
	double F = F_Res / (n + 1);
	double M_t = F * delta_XY;
	MatrixXd F_e(8, 1);
	int i;



	F_e.fill(0);
	F_e(1,0) = M_t;
	F_e(2,0) = F * 0.5;
	F_e(3,0) = F * l / 12;
	F_e(5,0) = M_t;
	F_e(6,0) = F * 0.5;
	F_e(7,0) = -F * l / 12;


	//transformierte Elementkraftvektoren
	MatrixXd FT(8, n);
	FT.fill(0);
		for (i = 0; i < n; i++){
		FT.block(0, 1 * i, 8, 1) = TMG.block(0, 8 * i, 8, 8).adjoint()*F_e.block(0, 0, 8, 1);
		
	}
	
	//Systemkraftvektor
	MatrixXd FG(4 * (n + 1), 1);
	MatrixXd FG_e(4 * (n + 1), 1);
	FG.fill(0);
	FG_e.fill(0);
	for (i = 0; i < n; i++){
		FG_e.fill(0);
		FG_e.block(4 * i, 0, 8, 1) = FT.block(0, 1 * i, 8, 1);
		FG = FG + FG_e;
		
	}

	VectorXd F_sys(4 * (n + 1));
	F_sys = FG;
	return(F_sys);
}












