/*----------------------------------------------------*/
/*--Author: Harrys Kon (Charalambos Konstantinou)-----*/
/*--W: https://harrys.fyi/----------------------------*/
/*--E: konharrys@gmail.com----------------------------*/
/*----------------------------------------------------*/
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <complex>

using namespace std;
using std::string;
using std::vector;
using std::cout;
using std::endl;

using vd = vector<double>;
using matr = vector<vd>;


void zdatas(int num, vd zd[], vd out[], matr G, matr B, matr bbus)
{
    static int offset = 0;
    //cout << "\n=======================\noffset=" << offset << '\n';

    string name = string()+"files/input"+std::to_string(num)+".txt";
    ifstream in(name);
    if ( !in ) throw "cannot open "+name;

    for ( int i = 0; i < offset; i++ )
    {
        string line;
        getline(in, line);
    }

    int size = 0;
    vd v[6];
    while (1)
    {
        string line;
        getline(in, line);
        if ( line.empty() || !in )
        {
            if ( size ) break;
            //cout << "NEW CYCLE [" << line << "]\n";
            offset = 0;
            in.close();
            zdatas(num,zd,out,G,B,bbus);
            return;
        }

        double a, b, c, d, e, f;
        istringstream is(line);
        is >> a >> b >> c >> d >> e >> f;
        if ( !is ) throw "input file is corrupted";

        if ( a == 1 && size > 0 ) break;

        v[0].push_back(a);
        v[1].push_back(b);
        v[2].push_back(c);
        v[3].push_back(d);
        v[4].push_back(e);
        v[5].push_back(f);

        size++;
    }

    for( int i=0; i<6; i++ )
    {
        zd[i].clear();
        zd[i].insert(zd[i].end(), v[i].begin(), v[i].end());
    }

    offset += size;
}
