#include<iostream>
#include"Element.h"
#include"GlobaData.h"
#include"Node.h"
#include"Siatka.h"
#include<iomanip>
#include<cmath>
#include<fstream>
using namespace std;
const double eps = 1e-12;

double** jakobian(Element el, int npc);	
double detJ(double **J);
double **OdwrotnoœæJacobianu(double **Jacobian, double det);
void wypiszJakobian(double** J);
void H_local_P_local(Element el, double k = 1, double alfa = 1);
double *GP(Siatka st);
double** GC(Siatka St);
double** GH(Siatka St);
void MES(Siatka St ,double **H, double**C,double*P, double tau , double t0, double simulation_time);
double **ludist(int n, double ** a);
double *lusolve(int n, double ** A, double * b, double *X);
int main()
{


	///////////////////////////////////////////////////////
	///////////////////////////////////////////////////////
	/////////////////////////////////////////////////////// DANE SIATKI

	double B = 0.100;
	double H = 0.100;
	int ilosc_wezlow_H = 31;
	int ilosc_wezlow_B = 31;
	int ilosc_elementow = (ilosc_wezlow_B - 1)*(ilosc_wezlow_H - 1);
	GlobaData Gd = GlobaData(B, H, ilosc_wezlow_B, ilosc_wezlow_H);

	Siatka St = Siatka(Gd.nb, Gd.nh, Gd.B, Gd.H);
	St.Stworz();

	//St.Show_All_Ell();
	//cout << endl;
	//cout << endl;
	// St.Show_All_Nodes();
	//cout << endl;
	//cout << endl;
	//St.Show_Status();
	/////////////////////////////////////////////////////////////// Zmienne test Case 1 
	
	double initial_temperature = 100;
	double simulation_time = 100;
	double simulation_step = 1;
	double ambient_temperature = 1200;
	double alfa = 300;
	double specific_heat = 700;
	double conductivity = 25;
	double denisity = 7800;


	///////////////////////////////////////////////////////////////PROGRAM G£OWNY
	double **Jacobian;
	for (int j = 0; j < ilosc_elementow; j++)
	{	
		//cout << "Element #" << j + 1 << endl;
		for (int i = 0; i < 4; i++)
		{
			
			Jacobian = jakobian(St.el[j], i);
			double det = detJ(Jacobian);
			double **Odwr = OdwrotnoœæJacobianu(Jacobian, det);
			St.el[j].f_dN_dXY(Odwr, i);
			
		}
		H_local_P_local(St.el[j],conductivity,alfa);
		St.el[j].f_P(ambient_temperature,alfa);
		St.el[j].f_C(specific_heat, denisity);
	}


	MES(St,GH(St), GC(St), GP(St), simulation_step, initial_temperature,simulation_time);




	system("PAUSE");
	return 0;
}

double** jakobian(Element el,int npc) {

	double **J = new double *[2];
	for (int i = 0; i < 2; i++)
		J[i] = new double[2];

	J[0][0] = el.dX_dKsi[npc];
	J[0][1] = el.dY_dKsi[npc];
	J[1][0] = el.dX_dEta[npc];
	J[1][1] = el.dY_dEta[npc];

	return J;
}
double detJ(double **J)
{
	double det = (J[0][0] * J[1][1]) - (J[0][1] * J[1][0]);

	return det;
}
double **OdwrotnoœæJacobianu(double **Jacobian , double det)
{
	double **J = new double *[2];
	for (int i = 0; i < 2; i++)
		J[i] = new double[2];



	J[0][0] = (1.0 / det) * Jacobian[1][1];
	J[0][1] = (1.0 / det) * (-Jacobian[0][1]);
	J[1][0] = (1.0 / det) * (-Jacobian[1][0]);
	J[1][1] = (1.0 / det) * Jacobian[0][0];

	return J;

}
void wypiszJakobian(double** J) {

	cout << ("Jakobian:") << endl;
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			cout << J[i][j] << "\t";
		}
		cout << endl;
	}
}
void H_local_P_local(Element el, double k, double alfa)
{
	double **HX = new double*[4];
	double **HY = new double*[4];
	double **calka = new double*[4];
	double **tab = new double*[4];
	double Jacob;
	double J[2][2];
	double det;

	J[0][0] = el.dX_dKsi[0];
	J[0][1] = el.dY_dKsi[0];
	J[1][0] = el.dX_dEta[0];
	J[1][1] = el.dY_dEta[0];

	det = (J[0][0] * J[1][1]) - (J[1][0] * J[0][1]);
	////////////////////////////wypelnianie zerami tabeli
	
	for (int i = 0; i < 4; i++)
	{
		HX[i] = new double[4];
		HY[i] = new double[4];
		calka[i] = new double[4];
		tab[i] = new double[4];
		for (int j = 0; j < 4; j++)
		{
			calka[i][j] = 0.0;
			el.H_P[j] = 0.0;
			el.H_local[i][j] = 0.0;
		}
	
	}
	///////////////////////////////////////////////////////////////// H_LOCAL
	for (int npc = 0; npc < 4; npc++)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				HX[j][i] = el.dNx_Npc[j][npc] * el.dNx_Npc[i][npc];
				HY[j][i] = el.dNy_Npc[j][npc] * el.dNy_Npc[i][npc];
			}
		}
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				el.H_local[i][j] += HX[i][j] + HY[i][j];
			}
		}
	}
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////H_P  "druga calka w równaniu"
	
	for (int i = 0; i < 4; i++)
	{
		
		if (el.Pow[i] == 1)
		{
			switch (i)
			{
			case 0:
				Jacob = (abs((el.nodes[0]->y - el.nodes[3]->y)) / 2.0);
				break;
			case 1:
				Jacob = (abs((el.nodes[0]->x - el.nodes[1]->x)) / 2.0);
				break;
			case 2:
				Jacob = (abs((el.nodes[1]->y - el.nodes[2]->y)) / 2.0);
				break;
			case 3:
				Jacob = (abs((el.nodes[3]->x - el.nodes[2]->x)) / 2.0);
				break;
			}
			for (int npc = 0; npc < 2; npc++)
			{
				tab[npc][0] = el.N1(el.K_npc_x[i][npc], el.K_npc_y[i][npc]);
				tab[npc][1] = el.N2(el.K_npc_x[i][npc], el.K_npc_y[i][npc]);
				tab[npc][2] = el.N3(el.K_npc_x[i][npc], el.K_npc_y[i][npc]);
				tab[npc][3] = el.N4(el.K_npc_x[i][npc], el.K_npc_y[i][npc]);
				//cout << tab[npc][0]<<"\t" << tab[npc][1] << "\t" << tab[npc][2] << "\t" << tab[npc][3];
				//cout << endl;	
			}
			//cout << endl;

			for (int k = 0; k < 4; k++)
			{
				for (int j = 0; j < 4; j++)
				{
					calka[j][k] = calka[j][k] + ((tab[0][k] * tab[0][j])*Jacob);
					calka[j][k] = calka[j][k] + ((tab[1][k] * tab[1][j])*Jacob);
				}
			}
		}
	}
	/*for (int k = 0; k < 4; k++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout<<calka[k][j] << "\t";
		}
		cout << endl;
	}*/
	////////////////////////wypis H_P i H_P *alfa
	for (int j = 0; j < 4; j++)
	{
		el.H_P[j] = el.H_P[j] * alfa;
		/*cout << el.H_P[j];
		cout << endl;*/

	}

	/////////////////////////////////////////////H_local Tworzenie
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			el.H_local[i][j] = (el.H_local[i][j] * det*k) + (calka[i][j])*alfa;
			//cout << el.H_local[i][j] << "\t";
		}
	//cout << endl;
	}
	//cout << endl;

	delete HX;
	delete HY;
	delete calka;
	delete tab;
}
double* GP(Siatka St)
{
	double matrix_size = St.nB*St.nH;
	double *P_Global = new double[matrix_size];
	for (int j = 0; j < matrix_size; j++)
	{
		P_Global[j] = 0.0;
	}
	for (int k = 0; k < St.ne; k++)
	{
		for (int j = 0; j < 4; j++)
		{
			P_Global[(St.el[k].Id[j]) - 1] += St.el[k].P_local[j];
			//cout << St.el[k].P_local[j] << "\t" << endl;
		}
	}
	//for (int j = 0; j < matrix_size; j++)
	//{
	//	cout << P_Global[j] << "\t" << endl;
	//}
	return P_Global;
}
double** GC(Siatka St)
{
	double matrix_size = St.nB*St.nH;
	double** GC = new double*[matrix_size];
	for (int i = 0; i < matrix_size; i++)
	{
		GC[i] = new double[matrix_size];
		for (int j = 0; j < matrix_size; j++)
		{
			GC[i][j] = 0.0;
		}
	
	}
	
	for (int k = 0; k < St.ne; k++)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				GC[(St.el[k].Id[i]) - 1][(St.el[k].Id[j]) - 1] += St.el[k].C_local[i][j];
			}
		}
	}
	//for (int i = 0; i < matrix_size; i++)
	//{
	//	for (int j = 0; j < matrix_size; j++)
	//	{
	//		cout << setw(14)<< GC[i][j];
	//	}
	//	cout << endl;
	//	cout << endl;
	//	cout << endl;
	//}
	return GC;
}
double** GH(Siatka St)
{
	double matrix_size = St.nB*St.nH;
	double **GH = new double *[matrix_size];
	for (int i = 0; i < matrix_size; i++)
	{
		GH[i] = new double[matrix_size];
		for (int j = 0; j < matrix_size; j++)
		{
			GH[i][j] = 0.0;
		}
	}
	for (int k = 0; k < St.ne; k++)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				GH[(St.el[k].Id[i]) - 1][(St.el[k].Id[j]) - 1] += St.el[k].H_local[i][j];
			}
		}
	}

	//for (int i = 0; i < matrix_size; i++)
	//{
	//	for (int j = 0; j < matrix_size; j++)
	//	{
	//		cout << setw(14)<< GH[i][j];
	//	}
	//	cout << endl;
	//	cout << endl;
	//	cout << endl;
	//}
	return GH;
}

void MES(Siatka St,double ** H, double ** C, double *P, double tau, double t0 , double simulation_time)
{	
	double matrix_size = St.nB*St.nH;
	double **MatrixH = new double*[matrix_size];
	double *P_Vector = new double[matrix_size];
	double **MatrixC = new double*[matrix_size];
	double *Temp_poczatkowa = new double[matrix_size];
	double *TempX = new double[matrix_size];
	double TempMax = 0;
	double TempMin = 0;
	int number = 0;
	int timer = tau;
	bool equal = false;
	int precision = 1000;
	fstream plik;
	for (int i = 0; i < matrix_size; i++)
	{
		MatrixH[i] = new double[matrix_size];
		P_Vector[i] = 0.0;
		Temp_poczatkowa[i] = t0;
	}

	for (int i = 0; i < matrix_size; i++)
	{
		MatrixC[i] = new double[matrix_size];
		for (int j = 0; j < matrix_size; j++)
		{
			MatrixH[i][j] = H[i][j] + (C[i][j] / tau);
			MatrixC[i][j] = C[i][j] / tau;

		}
	}
	/////////////////////////////////////////////////////////////////////// KOLEJNE KROKI CZASOWE
	plik.open("wyniki_ziemniak.txt", ios::out);
	if (!plik.good())
	{
		cout << "B³ad otwarcia!" << endl;
	}
	else
	{
		cout << setw(7) << "Time[s]\t" << setw(7) << "\tMinTemp[s]\t" << setw(7) << "MaxTemp[s]\t" << endl;
		while (equal != true)//for (double m = tau; m <= simulation_time; m += tau)
		{
			
			for (int i = 0; i < matrix_size; i++)
			{
				for (int j = 0; j < matrix_size; j++)
				{
					P_Vector[i] = P_Vector[i] + (MatrixC[i][j] * Temp_poczatkowa[j]);
					//cout << P_Vector[i] << "\t";
				}
				P_Vector[i] += P[i];
				//cout << P_Vector[i] << "\t";
			}
			//cout << endl;


			TempX = lusolve(matrix_size, ludist(matrix_size, MatrixH), P_Vector, TempX);

			TempMin = TempX[2];  //////randomowy indeks do znalezienia Tempmin
			for (int k = 0; k < matrix_size; k++)
			{
				if (TempX[k] > TempMax)
					TempMax = TempX[k];

				if (TempX[k] < TempMin)
					TempMin = TempX[k];
			}

			for (int j = 0; j < matrix_size; j++)
			{
				Temp_poczatkowa[j] = TempX[j];
				P_Vector[j] = 0.0;
			}
			cout << setw(7) << timer << "\t" << setw(7) << "\t" << TempMin << "\t\t" << setw(7) << TempMax << endl;
			number++;

			/////////////////////////////////////////////////////zapis do pliku 


			plik << timer << endl;
			for (int i = 0; i < St.nH; i++)
			{
				for (int j = i; j < matrix_size; j = j + St.nB)
				{
					plik << TempX[j] << "\t";
				}
				plik << endl;
			}
			plik << endl;

			equal = false;
			/*for (int i = 0; i < matrix_size; i++)
			{
				for (int j = 0; j < matrix_size; j++)
				{
					
					if (TempX[i] != TempX[j])
					{
						if ((round(TempX[i]*precision)) != (round(TempX[j] * precision)))
						{
							equal = false;
							break;
						}
					}
				}
				if (equal == false)
				{
					break;
				}
			}*/
			if ((round(TempMin * precision)) == (round(TempMax * precision)))
				equal = true;
			timer +=tau;
		}
	}
	plik.close();
	
	cout << "Struktura ziemniaka jest juz gotowa do spozycia -- Rozpoczynam proces zarumieniania" << endl;
}

// Funkcja dokonuje rozk³adu LU macierzy A
//----------------------------------------
double **ludist(int n, double ** a)
{	
	double **A = new double*[n];
	for (int i = 0; i < n; i++)
	{
		A[i] = new double[n];
		for (int j = 0; j < n; j++)
		{
			A[i][j] = a[i][j];
		}
	}

	int i, j, k;

	for (k = 0; k < n - 1; k++)
	{
		if (fabs(A[k][k]) < eps) return false;

		for (i = k + 1; i < n; i++)
			A[i][k] /= A[k][k];

		for (i = k + 1; i < n; i++)
			for (j = k + 1; j < n; j++)
				A[i][j] -= A[i][k] * A[k][j];
	}
	return A;
}

// Funkcja wyznacza wektor X na podstawie A i B
//---------------------------------------------
double *lusolve(int n, double ** A, double * b, double *x)
{
	double *B = new double[n];
	double *X = new double[n];
	for (int i = 0; i < n; i++)
	{
		
		B[i] = b[i];
		X[i] = x[i];
	}
	int    i, j;
	double s;

	X[0] = B[0];

	for (i = 1; i < n; i++)
	{
		s = 0;

		for (j = 0; j < i; j++) s += A[i][j] * X[j];

		X[i] = B[i] - s;
	}

	if (fabs(A[n - 1][n - 1]) < eps) return false;

	X[n - 1] /= A[n - 1][n - 1];

	for (i = n - 2; i >= 0; i--)
	{
		s = 0;

		for (j = i + 1; j < n; j++) s += A[i][j] * X[j];

		if (fabs(A[i][i]) < eps) return false;

		X[i] = (X[i] - s) / A[i][i];
	}

	
		//for (int i = 0; i < n; i++)
		//{
		//	cout << X[i] << "\t";
		//}
		//cout << endl;
	
	return X;
}