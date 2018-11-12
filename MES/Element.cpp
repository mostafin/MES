#include "Element.h"
#include"Node.h"
#include<cmath>
#include<iostream>
using namespace std;

Element::Element()
{
	this->nodes = new Node*[4];
	
	this->Npc[0][0] = -0.577350269;
	this->Npc[0][1] = -0.577350269;
	this->Npc[1][0] =  0.577350269;
	this->Npc[1][1] = -0.577350269;
	this->Npc[2][0] =  0.577350269;
	this->Npc[2][1] =  0.577350269;
	this->Npc[3][0] = -0.577350269;
	this->Npc[3][1] =  0.577350269;

	this->H_local = new double*[4];
	for (int i = 0; i < 4; i++)
		H_local[i] = new double[4];


	this->K_npc_x[0][0] = -1.0;
	this->K_npc_y[0][0] = -0.577350269;
	this->K_npc_x[0][1] = -1.0;
	this->K_npc_y[0][1] = 0.577350269;
	this->K_npc_x[1][0] = -0.577350269;
	this->K_npc_y[1][0] = -1.0;
	this->K_npc_x[1][1] = 0.577350269;
	this->K_npc_y[1][1] = -1.0;
	this->K_npc_x[2][0] = 1.0;
	this->K_npc_y[2][0] = -0.577350269;
	this->K_npc_x[2][1] = 1.0;
	this->K_npc_y[2][1] = 0.577350269;
	this->K_npc_x[3][0] = -0.577350269;
	this->K_npc_y[3][0] = 1.0;
	this->K_npc_x[3][1] = 0.577350269;
	this->K_npc_y[3][1] = 1.0;



}


Element::~Element()
{
}

double Element::dN1dKsi(double ksi, double eta)
{
	return -0.25 * (1 - eta);
}

double Element::dN2dKsi(double ksi, double eta)
{
	return 0.25 * (1 - eta);
}

double Element::dN3dKsi(double ksi, double eta)
{
	return 0.25 * (1 + eta);
}

double Element::dN4dKsi(double ksi, double eta)
{
	return -0.25 * (1 + eta);
}

double Element::dN1dEta(double ksi, double eta)
{
	return -0.25 * (1 - ksi);
}

double Element::dN2dEta(double ksi, double eta)
{
	return -0.25 * (1 + ksi);
}

double Element::dN3dEta(double ksi, double eta)
{
	return 0.25 * (1 + ksi);
}

double Element::dN4dEta(double ksi, double eta)
{
	return 0.25 * (1 - ksi);
}

double Element::N1(double ksi, double eta)
{
	return(0.25*(1 - ksi)*(1 - eta));
}

double Element::N2(double ksi, double eta)
{
	return (0.25*(1 + ksi)*(1 - eta));
}

double Element::N3(double ksi, double eta)
{
	return (0.25*(1 + ksi)*(1 + eta));
}

double Element::N4(double ksi, double eta)
{
	return (0.25*(1 - ksi)*(1 + eta));
}

void Element::f_dN_dXY(double ** OdwrotnoscJ, int npc)
{
	for (int i = 0; i < 4; i++)
	{
		this->dNx_Npc[i][npc] = (OdwrotnoscJ[0][0] * this->dN_dKsi[npc][i]) + (OdwrotnoscJ[0][1] * this->dN_dEta[npc][i]);
		this->dNy_Npc[i][npc] = (OdwrotnoscJ[1][0] * this->dN_dKsi[npc][i]) + (OdwrotnoscJ[1][1] * this->dN_dEta[npc][i]);
	}	
}

void Element::f_Pow()
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = i+1; j < 4; j++)
		{
			if (this->nodes[i]->status == 1 && this->nodes[j]->status == 1)
			{
				if (i == 0 && j == 3)
					this->Pow[0] = true;


				if (i == 0 && j == 1)
					this->Pow[1] = true;


				if (i == 1 && j == 2)
					this->Pow[2] = true;


				if (i == 2 && j == 3)
					this->Pow[3] = true;


			}
		}
	}
	for (int k = 0; k < 4; k++)
	{
		if (this->Pow[k] != true)
			this->Pow[k] = false;
	}
}

void Element::Dane()
{
	for (int i = 0; i < 4; i++)
	{
		this->dN_dKsi[i][0] = dN1dKsi(Npc[i][0], Npc[i][1]);
		this->dN_dKsi[i][1] = dN2dKsi(Npc[i][0], Npc[i][1]);
		this->dN_dKsi[i][2] = dN3dKsi(Npc[i][0], Npc[i][1]);
		this->dN_dKsi[i][3] = dN4dKsi(Npc[i][0], Npc[i][1]);

		this->dN_dEta[i][0] = dN1dEta(Npc[i][0], Npc[i][1]);
		this->dN_dEta[i][1] = dN2dEta(Npc[i][0], Npc[i][1]);
		this->dN_dEta[i][2] = dN3dEta(Npc[i][0], Npc[i][1]);
		this->dN_dEta[i][3] = dN4dEta(Npc[i][0], Npc[i][1]);
	}	 
	for (int i = 0; i < 4; i++)
	{
		this->dX_dKsi[i] = (dN_dKsi[i][0] * this->nodes[0]->x) + (dN_dKsi[i][1] * this->nodes[1]->x) + (dN_dKsi[i][2] * this->nodes[2]->x) + (dN_dKsi[i][3] * this->nodes[3]->x);
		this->dX_dEta[i] = (dN_dEta[i][0] * this->nodes[0]->x) + (dN_dEta[i][1] * this->nodes[1]->x) + (dN_dEta[i][2] * this->nodes[2]->x) + (dN_dEta[i][3] * this->nodes[3]->x);
		this->dY_dKsi[i] = (dN_dKsi[i][0] * this->nodes[0]->y) + (dN_dKsi[i][1] * this->nodes[1]->y) + (dN_dKsi[i][2] * this->nodes[2]->y) + (dN_dKsi[i][3] * this->nodes[3]->y);
		this->dY_dEta[i] = (dN_dEta[i][0] * this->nodes[0]->y) + (dN_dEta[i][1] * this->nodes[1]->y) + (dN_dEta[i][2] * this->nodes[2]->y) + (dN_dEta[i][3] * this->nodes[3]->y);
	}
	f_Pow();
}

void Element::f_P(double temp , double alfa) 
{
	double tab[2][4];
	double **P_local_buff = new double*[4];
	double Jacob;
	bool flag;
	/////////////////////////zerowanie wektora

	for (int i = 0; i < 4; i++)
	{
		P_local_buff[i] = new double[4];
		this->P_local[i] = 0.0;
	}
	for (int i = 0; i < 4; i++)
	{
		flag = false;
		if (this->Pow[i] == 1)
		{
			switch (i)
			{
			case 0:
				Jacob = (abs((this->nodes[0]->y - this->nodes[3]->y)) / 2.0);
				break;
			case 1:
				Jacob = (abs((this->nodes[0]->x - this->nodes[1]->x)) / 2.0);
				break;
			case 2:
				Jacob = (abs((this->nodes[1]->y - this->nodes[2]->y)) / 2.0);
				break;
			case 3:
				Jacob = (abs((this->nodes[3]->x - this->nodes[2]->x)) / 2.0);
				break;
			}
			for (int npc = 0; npc < 2; npc++)
			{
				tab[npc][0] = N1(this->K_npc_x[i][npc], this->K_npc_y[i][npc]);
				tab[npc][1] = N2(this->K_npc_x[i][npc], this->K_npc_y[i][npc]);
				tab[npc][2] = N3(this->K_npc_x[i][npc], this->K_npc_y[i][npc]);
				tab[npc][3] = N4(this->K_npc_x[i][npc], this->K_npc_y[i][npc]);
				//cout << tab[npc][0]<<"\t" << tab[npc][1] << "\t" << tab[npc][2] << "\t" << tab[npc][3];
				//cout << endl;	
			}
			//cout << endl;

			for (int j = 0; j < 4; j++)
			{
				P_local_buff[i][j] = (tab[0][j] + tab[1][j])*Jacob;		
			}

			flag = true;
		}
		if (flag == false)
		{
			for (int j = 0; j < 4; j++)
			{
				P_local_buff[i][j] = 0.0;
			}
		}
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			//cout << P_local_buff[i][j] << "\t";
		}
		//cout << endl;
	}
	//cout << endl;

	//////////////////////////wypis P_local i P_local *alfa
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			this->P_local[i] += P_local_buff[j][i];
		}
		this->P_local[i] = this->P_local[i] *alfa * temp;
	   // cout << this->P_local[i] << "\t";
	}
	//cout << endl;

	delete P_local_buff;
}

void Element::f_C(double c, double ro)
{
	double **tab = new double*[4];
	double **C_local_buff= new double*[4];
	///////////////////////////////////////////zerowanie macierzy C_lokal
	for (int i = 0; i < 4; i++)
	{
		tab[i] = new double[4];
		C_local_buff[i] = new double[4];
		for (int j = 0; j < 4; j++)
		{
			this->C_local[i][j] = 0.0;
			C_local_buff[i][j] = 0.0;
			
		}
	}

	for (int npc = 0; npc < 4; npc++)
	{
		tab[npc][0] = N1(this->Npc[npc][0], this->Npc[npc][1]);
		tab[npc][1] = N2(this->Npc[npc][0], this->Npc[npc][1]);
		tab[npc][2] = N3(this->Npc[npc][0], this->Npc[npc][1]);
		tab[npc][3] = N4(this->Npc[npc][0], this->Npc[npc][1]);

		//cout << tab[npc][0]<<"\t" << tab[npc][1]<<"\t "<< tab[npc][2] << "\t" << tab[npc][3] << endl;
	}
	//cout << endl;

	double **J = new double *[2];
	for (int i = 0; i < 2; i++)
		J[i] = new double[2];

	double *det = new double[4];
	for (int npc = 0; npc < 4; npc++)
	{
		det[npc] = 0.0;
		J[0][0] = this->dX_dKsi[npc];
		J[0][1] = this->dY_dKsi[npc];
		J[1][0] = this->dX_dEta[npc];
		J[1][1] = this->dY_dEta[npc];

		det[npc] = (J[0][0] * J[1][1]) - (J[1][0] * J[0][1]);
	/*	cout << J[0][0] << "\t " << J[0][1] << endl;
		cout << J[1][0] << "\t " << J[1][1] << endl;
		cout << endl;*/

	
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				this->C_local[i][j] = this->C_local[i][j] + (tab[npc][i] * tab[npc][j] * det[npc]);  
			}
		}
	}

	
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			this->C_local[i][j] = (this->C_local[i][j])* c * ro;
			//cout << C_local[i][j] << "\t";
		}
		//cout << endl;

	}
	delete J;
	delete tab;
	delete C_local_buff;
	delete det;
}
