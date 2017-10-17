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


void zdatas(int num, vd zd[], vd out[], int bus, matr G, matr B, matr bbus, bool hitlflag)
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
            zdatas(num,zd,out,bus,G,B,bbus,hitlflag);
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

    void update(vd * v, int num, vd out[], int bus, matr G, matr B, matr bbus);
    //if (hitlflag) 
        update(v,num,out,bus,G,B,bbus);

    for( int i=0; i<6; i++ )
    {
        zd[i].clear();
        zd[i].insert(zd[i].end(), v[i].begin(), v[i].end());
    }

    offset += size;
}

void update(vd * v, int num, vd out[], int bus, matr G, matr B, matr bbus)
{
    const int Vo = 0;
    const int Do = 1;

    if( out[Vo].empty() ){
        return;
    } 

    double input = 0;
    {
        std::ifstream in("files/wls.in");
        if( !in )
        {
            cout<<"No wls input !\r";
            return;
        }
        in >> input;
    }

    const int Type = 1;
    const int Val = 2;
    const int From = 3;
    const int To = 4;

    // here update code
    double dum = 0;
    unsigned int i,j;
    for ( size_t k = 0; k < v[0].size(); k++ )
    {
        if ( v[Type][k] == 1 && ( v[From][k] == bus && v[To][k] == 0))
        {
            v[Val][k] = input;
        }
        // Real Power Injection 
        if ( v[Type][k] == 2 && ( v[From][k] == bus && v[To][k] == 0))
        {
            dum = 0;
            for (j=0; j<out[Vo].size(); j++) {
                dum+=out[Vo][j]*(G[bus-1][j]*cos(out[Do][bus-1] - out[Do][j]) + B[bus-1][j]*sin(out[Do][bus-1] - out[Do][j]));
            }
            v[Val][k] = dum*input;
        }
        // Reactive Power Injection 
        if ( v[Type][k] == 3 && (v[From][k] == bus && v[To][k] == 0))
        {
            dum = 0;
            for (j=0; j<out[Vo].size(); j++) {
                dum+=out[Vo][j]*(G[bus-1][j]*sin(out[Do][bus-1] - out[Do][j]) - B[bus-1][j]*cos(out[Do][bus-1] - out[Do][j]));
            }
            v[Val][k] = dum*input;
        }
        // Real Power Flow
        if ( v[Type][k] == 4 && (v[From][k] == bus || v[To][k] == bus))
        {
            if (v[From][k] == bus) {
                j = int(v[To][k]-1);
                v[Val][k] = input*input*G[bus-1][j] - input*out[Vo][j]*(G[bus-1][j]*cos(out[Do][bus-1] - out[Do][j])  + B[bus-1][j]*sin(out[Do][bus-1] - out[Do][j]));
            }
            else if (v[To][k] == bus) {
                i = int(v[From][k]-1);
                v[Val][k] = out[Vo][i]*out[Vo][i]*G[i][bus-1] - out[Vo][i]*input*(G[i][bus-1]*cos(out[Do][i] - out[Do][bus-1])  + B[i][bus-1]*sin(out[Do][i] - out[Do][bus-1]));
            }
        }
        // Reactive Power Flow  
        if ( v[Type][k] == 5 && (v[From][k] == bus || v[To][k] == bus))
        {
            if (v[From][k] == bus) {
                j = int(v[To][k]-1);
                v[Val][k] = -input*input*(B[bus-1][j] + bbus[bus-1][j]) + input*out[Vo][j]*(B[bus-1][j]*cos(out[Do][bus-1] - out[Do][j])  - B[bus-1][j]*sin(out[Do][bus-1] - out[Do][j]));
            }
            else if (v[To][k] == bus) {
                i = int(v[From][k]-1);
                v[Val][k] = -out[Vo][i]*out[Vo][i]*(B[i][bus-1] + bbus[i][bus-1]) + input*out[Vo][i]*(B[i][bus-1]*cos(out[Do][i] - out[Do][bus-1])  - B[i][bus-1]*sin(out[Do][i] - out[Do][bus-1]));
            }
        }
    }
}

