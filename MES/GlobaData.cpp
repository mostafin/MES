#include "GlobaData.h"



GlobaData::GlobaData(double B, double H, double nb, double nh)
{
	this->B = B/(nb-1.0);
	this->H = H/(nh-1.0);
	this->nb = nb;
	this->nh = nh;
}


GlobaData::~GlobaData()
{
}
