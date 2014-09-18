#include "stdafx.h"
#include "up_konstanten.h"
#include "up_ode.h"
#include "up_kraftvektor.h"
#include <Eigen/LU>
#include <Eigen/SVD>

MatrixXd up_ode(MatrixXd &Y, double &F_Res, M_struct &M1){


	VectorXd F_system(4 * (n + 1));
	
	
	F_system = up_kraftvektor(F_Res, M1);

	MatrixXd A(4 * (n + 1), 1);
	MatrixXd B(4 * (n + 1), 4 * (n + 1));
	A = -M1.GM_sys*Y.col(1) + M1.KM_sys*Y.col(0) + F_system;
	B = M1.MM_sys;


	MatrixXd dy_i(4 * (n + 1), 2);
	dy_i.fill(0);
	dy_i.col(0) = Y.col(1);
	dy_i.col(1) = B.lu().solve(A);

	return(dy_i);
}

