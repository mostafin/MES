#pragma once
#include"Node.h"
#include"Element.h"
class Siatka
{
public:

	int ne;
	int nn;
	int nH;
	int nB;
	double B;
	double H;

	Node*nl;
	Element *el;

	Siatka(int ne , int nn, double B, double H);
	void Stworz();
	void Show_All_Ell();
	void Show_All_Nodes();
	void Show_Status();
	void Show_Pow();
	~Siatka();
};

