#include "Siatka.h"
#include<iostream>
#include<iomanip>
#include<cmath>
using namespace std;


Siatka::Siatka(int nB , int nH , double B , double H)
{
	this->B = B;
	this->H = H;
	this->nH = nH;
	this->nB = nB;

	this->ne = (nB-1)*(nH-1);
	this->nn = nB*nH;
	this->nl = new Node[nn];
	this->el = new Element[ne];
}

void Siatka::Stworz()
{
////////////////////////////////////////////////////////////////////////////////////nodes
		double t_x = 0.0;
		double t_y = 0.0;
		int temp = 1;
		for (int i = 0; i < nn; i++)
		{
			
			this->nl[i].x = t_x;
			this->nl[i].y = t_y;


			t_y = t_y + H;
			if (temp == nH)
				{
					t_x = t_x + B;
					t_y = 0.0;

					temp = 0;
				}

			temp++;

			if (nl[i].x == 0 || nl[i].y == 0 || nl[i].y > H*nH-(H+0.0000001) || nl[i].x > B*nB-(B+0.0000001))
				this->nl[i].status = 1;
			else
				this->nl[i].status = 0;
 		}
//////////////////////////////////////////////////////////////////////////////elementy
	int number_el = 0;
	int buff = 0;
	
		for (int j = 0; j < ne ; j++)
		{
			
			if (buff == nH - 1)
			{
				number_el++;
				buff = 0;
			}
		
			for (int k = 0; k < 4; k++)
			{
				switch (k) {
				case 0:
					this->el[j].Id[0] = number_el + 1;
					this->el[j].nodes[0] = &nl[number_el];
					break;

				case 1:
					this->el[j].Id[1] = number_el + nH + 1;
					this->el[j].nodes[1] = &nl[number_el + nH];
					break;
				case 2:
					this->el[j].Id[2] = number_el + nH + 2;
					this->el[j].nodes[2] = &nl[number_el + nH + 1];
					break;

				case 3:
					this->el[j].Id[3] = number_el + 2;
					this->el[j].nodes[3] = &nl[number_el + 1];
					break;
				}
				
			}
			this->el[j].Dane();
			number_el++;
			buff++;
		}
		//////////////////////////////////////
	}

void Siatka::Show_All_Ell()
{

	for (int j = nH - 2; j >= 0; j--)
	{
		int  k;
		for (int i = j , k = 0; k < nB-1 ; i = i + (nB - 1), k++)
		{
			std::cout << std::setw(3) << this->el[i].Id[3] << "|" << std::setw(3) << this->el[i].Id[2] << std::setw(2) << "";
		}
		std::cout << std::endl;
		std::cout << "-------------------------------------------------------------------------------" << std::endl;
		for (int i = j , k = 0 ; k < nB-1 ; i = i + (nB - 1) , k++)
		{
			std::cout << std::setw(3) << this->el[i].Id[0] << "|" << std::setw(3) << this->el[i].Id[1] << std::setw(2) << "";
		}
		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
	}
		
}

void Siatka::Show_All_Nodes()
{
	for (int j = nH - 1; j >= 0; j--)
	{
		int  k;
		for (int i = j, k = 0; k < nB ; i = i + nB, k++)
		{
			std::cout << std::setw(3) <<"[" << this->nl[i].x << "," << std::setw(3) << this->nl[i].y << "] | ";
			//std::cout << std::setw(3) << "[" << i+1 << "] | ";
		}
		std::cout << std::endl;
		std::cout << "------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
		
	}
}

void Siatka::Show_Status()
{
	for (int j = nH - 1; j >= 0; j--)
	{
		int  k;
		for (int i = j, k = 0; k < nB; i = i + nB, k++)
		{
			std::cout << std::setw(2) << this->nl[i].status << "|";
		}
		std::cout << std::endl;
		std::cout << "------------------------------" << std::endl;
	}
}

void Siatka::Show_Pow()
{
	
}


Siatka::~Siatka()
{
}
