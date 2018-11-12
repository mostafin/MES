#pragma once
class GlobaData
{
public:

	double B;	//szerokosc
	double H;	//wysokosc
	double nb;	//ilosc wezlow po szerokosci
	double nh;	//ilosc wezlow po wysokosci
	GlobaData(double B, double H, double nb, double nh);
	~GlobaData();
};

