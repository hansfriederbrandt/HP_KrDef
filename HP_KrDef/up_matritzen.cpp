#include "stdafx.h"
#include "up_konstanten.h"
#include "up_matritzen.h"


void up_matritzen(M_struct &M1){


	MatrixXd KM(8, 8);
	MatrixXd KMBieg(4, 4);
	MatrixXd KMTors(2, 2);
	MatrixXd KMZug(2, 2);
	MatrixXd MM(8, 8);


	double Tors_Fakt;
	double ZG_Fakt;
	double Bieg_Fakt;
	double Massen_Fakt;
	double phi;
	int i;

	//Elementsteifigkeitsmatrix
	Bieg_Fakt = (EJ / pow(l, 3.0));
	KMBieg << 12, 6 * l, -12, 6 * l,
		6 * l, 4 * pow(l, 2.0), -6 * l, 2 * pow(l, 2.0),
		-12, -6 * l, 12, -6 * l,
		6 * l, 2 * pow(l, 2.0), -6 * l, 4 * pow(l, 2.0);
	KMBieg = KMBieg * Bieg_Fakt;

	//Torsionsmatrix
	Tors_Fakt = (G*It / l);
	KMTors << 1, -1,
		-1, 1;
	KMTors = KMTors *Tors_Fakt;

	//Elementsteifigkeitsmatrix Zug / Druck
	ZG_Fakt = (EA / l);
	KMZug << 1, -1,
		-1, 1;
	KMZug = KMZug*ZG_Fakt;

	//Elementsteifigkeitsmatrix für ein Element erstellen
	KM.fill(0);

	KM.block(2, 2, 2, 2) = KMBieg.block(0, 0, 2, 2);
	KM.block(6, 2, 2, 2) = KMBieg.block(2, 0, 2, 2);
	KM.block(2, 6, 2, 2) = KMBieg.block(0, 2, 2, 2);
	KM.block(6, 6, 2, 2) = KMBieg.block(2, 2, 2, 2);

	KM(0, 0) = KMZug(0, 0);
	KM(4, 0) = KMZug(1, 0);
	KM(0, 4) = KMZug(0, 1);
	KM(4, 4) = KMZug(1, 1);

	KM(1, 1) = KMTors(0, 0);
	KM(5, 1) = KMTors(1, 0);
	KM(1, 5) = KMTors(0, 1);
	KM(5, 5) = KMTors(1, 1);


	//Elementmassenmatrix

	MM << 140, 0, 0, 0, 70, 0, 0, 0,
		0, 140 * It, 0, 0, 0, 70 * It, 0, 0,
		0, 0, 156, -22 * l, 0, 0, 54, 13 * l,
		0, 0, -22 * l, 4 * pow(l, 2.0), 0, 0, -13 * l, -3 * pow(l, 2.0),
		70, 0, 0, 0, 140, 0, 0, 0,
		0, 70 * It, 0, 0, 0, 140 * It, 0, 0,
		0, 0, 54, -13 * l, 0, 0, 156, 22 * l,
		0, 0, 13 * l, -3 * pow(l, 2.0), 0, 0, 22 * l, 4 * pow(l, 2.0);

	Massen_Fakt = (rhoA*l / 420);
	MM = MM * Massen_Fakt;

	MatrixXd TM(8, 8);
	MatrixXd TMG(8, 8 * n);
	TM.fill(0);
	TMG.fill(0);

	MatrixXd TMtest(8, 16);
	TMtest.fill(0);

	for (i = 0; i < n; i++){





		phi = (i)*(2 - ss)*pi / n;

		TM << cos(phi), 0, sin(phi), 0, 0, 0, 0, 0,
			0, 1, 0, 0, 0, 0, 0, 0,
			-sin(phi), 0, cos(phi), 0, 0, 0, 0, 0,
			0, 0, 0, 1, 0, 0, 0, 0,
			0, 0, 0, 0, cos(phi), 0, sin(phi), 0,
			0, 0, 0, 0, 0, 1, 0, 0,
			0, 0, 0, 0, -sin(phi), 0, cos(phi), 0,
			0, 0, 0, 0, 0, 0, 0, 1;

		TMG.block(0, 8 * (i), 8, 8) = TM;

	}




	//Transformation der Elementsteifigkeitsmatrix
	MatrixXd KMT(8, 8 * n);
	KMT.fill(0);
	for (i = 0; i < n; i++){
		KMT.block(0, 8 * i, 8, 8) = TMG.block(0, 8 * i, 8, 8).adjoint() * KM * TMG.block(0, 8 * i, 8, 8);
	}


	//Transformation der Elementmassenmatrix
	MatrixXd MMT(8, 8 * n);
	MMT.fill(0);
	for (i = 0; i < n; i++){
		MMT.block(0, 8 * i, 8, 8) = TMG.block(0, 8 * i, 8, 8).adjoint() * MM * TMG.block(0, 8 * i, 8, 8);
	}


	//Systemsteifigkeitsmatrix
	MatrixXd KMG(4 * (n + 1), 4 * (n + 1));
	MatrixXd KMG_e(4 * (n + 1), 4 * (n + 1));
	KMG.fill(0);
	KMG_e.fill(0);
	for (i = 0; i < n; i++){
		KMG_e.fill(0);
		KMG_e.block(4 * i, 4 * i, 8, 8) = KMT.block(0, 8 * i, 8, 8);
		KMG += KMG_e;
	}

	//Systemmassenmatrix
	MatrixXd MMG(4 * (n + 1), 4 * (n + 1));
	MatrixXd MMG_e(4 * (n + 1), 4 * (n + 1));
	MMG.fill(0);
	MMG_e.fill(0);
	for (i = 0; i < n; i++){
		MMG_e.fill(0);
		MMG_e.block(4 * i, 4 * i, 8, 8) = MMT.block(0, 8 * i, 8, 8);
		MMG += MMG_e;
	}

	//Gesamtdämpfungsmatrix
	MatrixXd GMG(4 * (n + 1), 4 * (n + 1));
	GMG.fill(0);
	double alpha = 0;
	double beta = 0;
	GMG = alpha*MMG + beta*KMG;      //Dämpfungsmatrix Gamma(in Ns / m)


	M1.KM_sys = KMG;
	M1.MM_sys = MMG;
	M1.TM_sys = TMG;
	M1.GM_sys = GMG;

}

