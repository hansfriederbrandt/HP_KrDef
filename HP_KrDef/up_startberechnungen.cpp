#include "stdafx.h"
#include "up_konstanten.h"
#include "up_startberechnungen.h"

void up_startberechnungen(Vector2d &PS)
{
	
	//lokale Variable initialisieren
	MatrixXd PRand(6, 2);
	PRand << -0.037562690734863, 0.15856578063965,
		-0.037143062591553, 0.15856578063965,
		-0.034734931945801, 0.15868771362305,
		-0.034734931945801, 0.16132368469238,
		-0.037154541015625, 0.16142205810547,
		-0.0374965400695, 0.16142205810547;

	PRand = PRand * 1000;

	
	MatrixXd PS_e(P_Anz, 2);
	PS_e(0, 0) = 0, 0;

	MatrixXd Ixxyy_e(P_Anz, 2);
	Ixxyy_e(0, 0) = 0, 0;

	Vector2d Ixxyy;
	Ixxyy << 0, 0;

	double Ixy;

	Vector2d P1;
	P1 << 0, 0;

	Vector2d P2;
	P2 << 0, 0;

	VectorXd A_e(P_Anz);
	A_e(0) = 0;
	VectorXd  Ixy_e(P_Anz);


	double test[6] = { 0, 0, 0, 0, 0, 0 };
	double test2[6] = { 0, 0, 0, 0, 0, 0 };

	int j, i;



	//Koordinate Prand(k, 0) = x - Koordinate
	//Koordinate Prand(k, 1) = y - Koordinate
	for (i = 0; i < P_Anz; i++){

		if (i < (P_Anz - 1))
			j = i + 1;
		else
			j = 0;

		A_e(i) = (PRand(i, 0) * PRand(j, 1) - PRand(j, 0) * PRand(i, 1));

		PS_e(i, 0) = (PRand(i, 0) + PRand(j, 0)) * (PRand(i, 0) * PRand(j, 1) - PRand(j, 0) * PRand(i, 1));
		PS_e(i, 1) = (PRand(i, 1) + PRand(j, 1)) * (PRand(i, 0) * PRand(j, 1) - PRand(j, 0) * PRand(i, 1));

		//A = (A + A_e);

		/*PS_e_sum(0) = PS_e_sum(0) + PS_e(0);
		PS_e_sum(1) = PS_e_sum(1) + PS_e(1);*/

	}

	A = A_e.sum() / 2;
	PS = PS_e.leftCols(2).colwise().sum() / (6 * A);

	//Bestimmung des Flächenträgheitselemente unter Anwendung des Green'schen
	//Theorems bzw.Satz von Steiner

	for (i = 0; i < P_Anz; i++){
		if (i < (P_Anz - 1))
			j = i + 1;
		else
			j = 0;


		P1(0) = PRand(i, 0) - PS(0);
		P1(1) = PRand(i, 1) - PS(1);

		P2(0) = PRand(j, 0) - PS(0);
		P2(1) = PRand(j, 1) - PS(1);



		Ixxyy_e(i, 0) = (P2(0) - P1(0)) * (pow(P1(1), 3.0) + pow(P1(1), 2.0) * P2(1) + P1(1) * pow(P2(1), 2.0) + pow(P2(1), 3.0));
		Ixxyy_e(i, 1) = (P2(1) - P1(1)) * (pow(P1(0), 3.0) + pow(P1(0), 2.0) * P2(0) + P1(0) * pow(P2(0), 2.0) + pow(P2(0), 3.0));
		Ixy_e(i) = (P1(1) * P2(1) + 2 * P1(0) * P1(1) + 2 * P2(0) * P2(1) + P2(0) + P1(1)) * (P1(0) * P2(1) - P2(0) * P1(1));

	}

	Ixxyy = Ixxyy_e.colwise().sum() / (-12);
	Ixy = abs(Ixy_e.sum() / (-24));

	Ixxyy(0) = abs(Ixxyy(0));
	Ixxyy(1) = abs(Ixxyy(1));


	Jx = Ixxyy(0);     //[mm ^ 4] Querschnitts - Flächenträgheitsmoment in X - Richtung
	Jy = Ixxyy(1);     //[mm ^ 4] Querschnitts - Flächenträgheitsmoment in Y - Richtung

	Vector2d b_h_max;
	b_h_max(0) = PRand.col(0).maxCoeff() - PRand.col(0).minCoeff();
	b_h_max(1) = PRand.col(1).maxCoeff() - PRand.col(1).minCoeff();

	//%Abstände der prohjezierten Randflächen Mittelpunkte zum Schwerkraftmittelpunkt
	 delta_XY = PS(1) - ((PRand.col(1).maxCoeff() - PRand.col(1).minCoeff()) / 2 + PRand.col(1).minCoeff());

	double l_Balk = 5 * b_h_max.maxCoeff(); //[mm] Länge der einzelnen Balkenelemente für bestes Verhältnis für FEM    //n = input('Balkenanzahl n = ');   //Anzahl der Balkenelemente in die Ring eingeteilt wird.
	n = L / l_Balk; //[]; Anzahl der Balkenelemente
	l = L / n;          //[mm]; Elementlänge

	EJ = E*Jy;          //[Nmm ^ 2]
	EA = E*A;           //[N]
	rhoA = rho*A;     //[kg / mm]

}