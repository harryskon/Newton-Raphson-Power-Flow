/*----------------------------------------------------*/
/*--Author: Harrys Kon (Charalambos Konstantinou)-----*/
/*--W: https://harrys.fyi/----------------------------*/
/*--E: konharrys@gmail.com----------------------------*/
/*----------------------------------------------------*/
struct SysTable
{
	enum Typs { S14, S30 } typs; // system
	enum Typt { BUS, LINE } typt; // table
	const double * data;
	int sz, cnum;

	SysTable(Typs,Typt);
	int size() const { return sz/cnum; }
	int cols() const { return cnum; }
	double operator()(int c, int r) const { return data[(c-1)+cnum*(r-1)]; }
};

struct SysData
{
	SysTable bus, line;
	SysData(SysTable::Typs t): bus(t,SysTable::BUS), line(t,SysTable::LINE) {}
};