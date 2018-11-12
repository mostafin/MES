#pragma once
#include"Node.h"
class Element
{
public:

	int Id[4];
	bool Pow[4];
	double dN_dKsi[4][4];
	double dN_dEta[4][4];
	double Npc[4][2];
	double K_npc_x[4][2];
	double K_npc_y[4][2];
	double **H_local;
	double dX_dKsi[4];
	double dX_dEta[4];
	double dY_dKsi[4];
	double dY_dEta[4];
	double dNx_Npc[4][4];
	double dNy_Npc[4][4];
	double P_local[4];
	double C_local[4][4];
	double H_P[4];
	Node **nodes;

	Element();
	~Element();
	double dN1dKsi(double ksi, double eta);
	double dN2dKsi(double ksi, double eta);
	double dN3dKsi(double ksi, double eta);
	double dN4dKsi(double ksi, double eta);
	double dN1dEta(double ksi, double eta);
	double dN2dEta(double ksi, double eta);
	double dN3dEta(double ksi, double eta);
	double dN4dEta(double ksi, double eta);
	double N1(double ksi, double eta);
	double N2(double ksi, double eta);
	double N3(double ksi, double eta);
	double N4(double ksi, double eta);
	void f_dN_dXY(double **OdwrotnoscJ, int node);
	void f_Pow();
	void Dane();
	void f_P(double temp = 1, double alfa = 1);
	void f_C(double c = 1, double ro = 1);
};

