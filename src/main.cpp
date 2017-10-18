/*----------------------------------------------------*/
/*--Author: Harrys Kon (Charalambos Konstantinou)-----*/
/*--W: https://harrys.fyi/----------------------------*/
/*--E: konharrys@gmail.com----------------------------*/
/*----------------------------------------------------*/
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>
#include <functional>
#include <iomanip>

#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include "Matrix.h"
#include "systab.h"
 
using std::cout;
using namespace std::complex_literals;
using std::vector;
using std::cerr;
using std::endl;
using std::string;
using std::complex;
using std::transform;
using std::plus;
using std::minus;
using vd = vector<double>;
using vc = vector<complex<double>>;
using matr = vector<vd>;

void cinv( Matrix RealA, Matrix ImagA, Matrix& RealAinv, Matrix& ImagAinv);

void subtract(const vd &a, const vd &b, vd &c)
{
    transform(a.begin(),a.end(),b.begin(),c.begin(), [] (double a, double b) { return a - b; });
}

void add(const vd &a, const vd &b, vd &c)
{
    transform(a.begin(),a.end(),b.begin(),c.begin(), [] (double a, double b) { return a + b; });
}

void multiply(const vd &a, const vd &b, vd &c)
{
    transform(a.begin(),a.end(),b.begin(),c.begin(), [] (double a, double b) { return a*b; });
}

double matrmultiply(matr A, matr B, int i, int j)
{
    double val= 0;
    for (unsigned int k = 0; k < A[0].size(); k++)
        val += (A[i][k])*(B[k][j]);
    return val;
}

void vecmultiply(matr A, vd &B, vd &C)
{
    C.clear();
    double val;
    for (unsigned int i = 0; i < A.size(); i++){
        val = 0;
        for (unsigned int j = 0; j < A[0].size(); j++)
            val += (A[i][j])*(B[j]);
        C.push_back(val);
    }
}

void cvecmultiply(vector<vc> A, vc &B, vc &C)
{
    C.clear();
    std::complex<double> val;
    for (unsigned int i = 0; i < A.size(); i++){
        val = 0. + 0i;
        for (unsigned int j = 0; j < A[0].size(); j++)
            val += (A[i][j])*(B[j]);
        C.push_back(val);
    }
}

void printMatrix(matr A) 
{
    for (unsigned int i=0; i<A.size(); i++) {
        for (unsigned int j=0; j<A[0].size(); j++)
            cout << A[i][j] << "\t";
        cout << "\n";
    }
    cout << endl;
}

void cprintMatrix(vector<vc> A) 
{
    for (unsigned int i=0; i<A.size(); i++) {
        for (unsigned int j=0; j<A[0].size(); j++)
            cout << A[i][j] << "\t";
        cout << "\n";
    }
    cout << endl;
}

void printVector(vd A) 
{
    for (unsigned int j=0; j<A.size(); j++)
        cout << A[j] << "\n";
    cout << endl;
}

void cprintVector(vc A) 
{
    for (unsigned int j=0; j<A.size(); j++)
        cout << A[j] << "\n";
    cout << endl;
}

double sumMatrix(matr A)
{
    double s=0;
    for (unsigned int i=0; i<A.size(); i++) {
        for (unsigned int j=0; j<A[0].size(); j++)
            s+=A[i][j];
    }
    return s;    
}

double sumVector(vd A)
{
    double s=0;
    for (unsigned int j=0; j<A.size(); j++)
        s+=A[j];
    return s;    
}

std::complex<double> csumMatrix(vector<vc> A)
{
    std::complex<double> s = 0. + 0i;
    for (unsigned int i=0; i<A.size(); i++) {
        for (unsigned int j=0; j<A[0].size(); j++)
            s+=A[i][j];
    }
    return s;    
}

std::complex<double> csumVector(vc A)
{
    std::complex<double> s = 0. + 0i;
    for (unsigned int j=0; j<A.size(); j++)
        s+=A[j];
    return s;    
}


int main(int argc, char* argv[])
try 
{
    // Check the number of parameters
    if (argc == 1) {
        cout << "Newton Raphson Power Flow of IEEE " << argv[1] << "-bus system" << endl;
    } 
    else if (argc < 1) {
        cerr << "Usage: " << argv[0] << " <number of IEEE bus system>" << endl;
        return 1;
    }
    
    
    SysData s14(SysTable::S14);
    SysData s30(SysTable::S30);
    SysData * sn;

    unsigned int N=atoi(argv[1]);
    const unsigned int bMVA = 100;
    const std::complex<double> cbMVA = bMVA;

    if (N==14)
        sn = &s14;
    else if (N==30)
        sn = &s30;
    else
        throw "Bad input IEEE system argument";

    
    //Form Admittance (Y) And Impedance (Z) (inv(ybus)) Bus Formation
    unsigned int i,j,m,n;
    unsigned int size_z = sn->line.size();
    vd r,x,fb,tb,a;
    vc z,y,b;
    matr bbus(N,vd(N,0)), G(N,vd(N)), B(N,vd(N));
    vector<vc> ybus(N,vc(N,0));

    for(i=0;i<size_z;i++) {
        fb.push_back(sn->line(1,i+1));
        tb.push_back(sn->line(2,i+1));
        r.push_back(sn->line(3,i+1));
        x.push_back(sn->line(4,i+1));
        b.push_back((sn->line(5,i+1))*1i);
        a.push_back(sn->line(6,i+1));
        z.push_back(r[i] + x[i]*1i);
        y.push_back(1./z[i]);  
    }

    //Formation of the off Diagonal Elements of Y-bus and B-bus
    for(i=0;i<size_z;i++) {
        ybus[int(fb[i])-1][int(tb[i])-1] = ybus[int(fb[i])-1][int(tb[i])-1] - y[i]/a[i];
        ybus[int(tb[i])-1][int(fb[i])-1] = ybus[int(fb[i])-1][int(tb[i])-1];

        bbus[int(fb[i])-1][int(tb[i])-1] = b[i].imag();
        bbus[int(tb[i])-1][int(fb[i])-1] = bbus[int(fb[i])-1][int(tb[i])-1];
    }

    // Formation of Diagonal Elements of Y-bus
    for(i=0;i<N;i++) {
        for(j=0;j<size_z;j++){
            if ((unsigned)((fb[j]) - 1) == i) {
                ybus[i][i] = ybus[i][i] + y[j]/(a[j]*a[j]) + b[j];
            }
            else if ((unsigned)((tb[j]) - 1) == i) {
                ybus[i][i] = ybus[i][i] + y[j] + b[j];
            }
        }
    }

    for(i=0;i<N;i++) {
        for(j=0;j<N;j++){
            G[i][j] = ybus[i][j].real();
            B[i][j] = ybus[i][j].imag();    
        }
    }

    
    // Get bus data
    unsigned int size_b = sn->bus.size();
    vd bus, typeb, V, del, Pg, Qg, Pl, Ql, Qmin, Qmax;

    for(i=0;i<size_b;i++) {
        bus.push_back(sn->bus(1,i+1));
        typeb.push_back(sn->bus(2,i+1));
        V.push_back(sn->bus(3,i+1));
        del.push_back(sn->bus(4,i+1));
        Pg.push_back((sn->bus(5,i+1))/bMVA);
        Qg.push_back((sn->bus(6,i+1))/bMVA);
        Pl.push_back((sn->bus(7,i+1))/bMVA);
        Ql.push_back((sn->bus(8,i+1))/bMVA);
        Qmin.push_back((sn->bus(9,i+1))/bMVA);
        Qmax.push_back((sn->bus(10,i+1))/bMVA); 
    }

    vd Psp(Pg.size());
    vd Qsp(Qg.size());

    subtract(Pg,Pl,Psp);
    subtract(Qg,Ql,Qsp);


    vd pv,pq;

    for(i=0; i<typeb.size(); i++) {
        if ((typeb[i] == 1) || (typeb[i] == 2))
            pv.push_back(bus[i]);
        else if (typeb[i] == 3)
            pq.push_back(bus[i]);
        else
            throw "bus file is corrupted";
    }

    unsigned int iter = 1;
    double tol = 5.0;
    double QG;

    while (tol > 1e-5) {
        
        vd P(N,0), Q(N,0);

        // Calculate P and Q
        for(i=0;i<N;i++) {
            for(j=0;j<N;j++){
                P[i] = P[i] + V[i]*V[j]*(G[i][j]*cos(del[i]-del[j]) + B[i][j]*sin(del[i]-del[j]));
                Q[i] = Q[i] + V[i]*V[j]*(G[i][j]*sin(del[i]-del[j]) - B[i][j]*cos(del[i]-del[j]));
            }
        }

        // Check Q-limit violations
        if ((iter > 2) && (iter <= 7)) {
            for(i=1;i<N;i++) {
                if (typeb[i] == 2) {
                    QG = Q[i] + Ql[i];
                    if (QG < Qmin[i])
                        V[i]+=0.01;
                    else if (QG > Qmax[i])
                        V[i]-=0.01;
                }
            }
        }

        // Calculate change from specified value
        vd dPa(P.size());
        vd dQa(Q.size());

        subtract(Psp,P,dPa);
        subtract(Qsp,Q,dQa);

        unsigned int k=0;
        vd dQ(pq.size(),0);
        vd dP;

        for (i=0; i<N; i++) {
            if (typeb[i] == 3) {
                dQ[k] = dQa[i];
                k+=1;
            }
        }
 
        dP.insert(dP.end(), dPa.begin()+1, dPa.end());

        // Mismatch vector
        vd M;

        M.insert(M.end(), dP.begin(), dP.end());
        M.insert(M.end(), dQ.begin(), dQ.end());

        // Jacobian..
        
        matr J1(N-1,vd(N-1,0));
        matr J2(N-1,vd(pq.size(),0));
        matr J3(pq.size(),vd(N-1,0));
        matr J4(pq.size(),vd(pq.size(),0));


        i = J1.size()+J3.size();            
        matr J(i);
        
       // J1 - Derivative of Real Power Injection with respect to angles.
        for(i=0;i<N-1;i++) {
            m = i+1;
            for(j=0;j<N-1;j++){
                n = j+1;
                if (m == n) {
                    for(k=0;k<N;k++){
                        J1[i][j] = J1[i][j] + V[m]*V[k]*(-G[m][k]*sin(del[m]-del[k]) + B[m][k]*cos(del[m]-del[k]));
                    }
                    J1[i][j] = J1[i][j] - V[m]*V[m]*B[m][m];  
                }
                else {
                    J1[i][j] = V[m]*V[n]*(G[m][n]*sin(del[m]-del[n]) - B[m][n]*cos(del[m]-del[n]));
                }    
            }
        }

        for(i=0;i<N-1;i++) {
            m = i+1;
            for(j=0;j<pq.size();j++){
                n = pq[j] -1;
                if (n < 0) n = 0;
                if (m == n) {
                    for(k=0;k<N;k++){
                        J2[i][j] = J2[i][j] + V[k]*(G[m][k]*cos(del[m]-del[k]) + B[m][k]*sin(del[m]-del[k]));
                    }
                    J2[i][j] = J2[i][j] + V[m]*G[m][m];  
                }
                else {
                    J2[i][j] = V[m]*(G[m][n]*cos(del[m]-del[n]) + B[m][n]*sin(del[m]-del[n]));
                }    
            }
        }

        for(i=0;i<pq.size();i++) {
            m = pq[i] - 1;
            if (m < 0) m = 0;
            for(j=0;j<N-1;j++){
                n = j+1;
                if (m == n) {
                    for(k=0;k<N;k++){
                        J3[i][j] = J3[i][j] + V[m]*V[k]*(G[m][k]*cos(del[m]-del[k]) + B[m][k]*sin(del[m]-del[k]));
                    }
                    J3[i][j] = J3[i][j] - V[m]*V[m]*G[m][m];  
                }
                else {
                    J3[i][j] = V[m]*V[n]*(-G[m][n]*cos(del[m]-del[n]) - B[m][n]*sin(del[m]-del[n]));
                }    
            }
        }

        for(i=0;i<pq.size();i++) {
            m = pq[i] - 1;
            if (m < 0) m = 0;
            for(j=0;j<pq.size();j++){
                n = pq[j] -1;
                if (n < 0) n = 0;
                if (m == n) {
                    for(k=0;k<N;k++){
                        J4[i][j] = J4[i][j] + V[k]*(G[m][k]*sin(del[m]-del[k]) - B[m][k]*cos(del[m]-del[k]));
                    }
                    J4[i][j] = J4[i][j] - V[m]*B[m][m];  
                }
                else {
                    J4[i][j] = V[m]*(G[m][n]*sin(del[m]-del[n]) - B[m][n]*cos(del[m]-del[n]));
                }    
            }
        }

        // Jacobian J
        for(i=0;i<J1.size();i++) {
            J[i].insert(J[i].end(), J1[i].begin(), J1[i].end());
            J[i].insert(J[i].end(), J2[i].begin(), J2[i].end());
        }
        unsigned int tmp = J1.size();
        for(i=0;i<J3.size();i++) {
            J[tmp+i].insert(J[tmp+i].end(), J3[i].begin(), J3[i].end());
            J[tmp+i].insert(J[tmp+i].end(), J4[i].begin(), J4[i].end());
        }
        
        // calc inv(J)
        unsigned int size_J = J.size();

        Matrix RealA1(size_J,size_J);
        Matrix ImagA1(size_J,size_J);
        Matrix RealA1inv(size_J,size_J);
        Matrix ImagA1inv(size_J,size_J);

        for (i=0; i<size_J; i++){
            for (j=0; j<size_J; j++){
                RealA1(i+1,j+1) = J[i][j];
                ImagA1(i+1,j+1) = 0;
            }   
        }
        cinv(RealA1, ImagA1, RealA1inv, ImagA1inv);

        matr invJ(size_J,vd(size_J));

        for (i=0; i<size_J; i++){
            for (j=0; j<size_J; j++){
                invJ[i][j] = RealA1inv(i+1,j+1);
            }   
        }

        //Correction vector
        vd X;
        vecmultiply(invJ,M,X);

        // Change in voltage angle & voltage magnitude
        vd dTh,dV;
        dTh.insert(dTh.end(), X.begin(), X.begin()+N-1);
        dV.insert(dV.end(), X.begin()+N-1, X.end());

        // Updating state vectors

        // Voltage angle
        for (i=1; i<N; i++){
            del[i]+=dTh[i-1];
        }

        // Voltage magnitude
        k = 0;
        for (i=1; i<N; i++) {
            if (typeb[i] == 3) {
                V[i]+=dV[k];
                k+=1;
            }
        }

        vd absM(M.size());
        for (i = 0; i < M.size(); i++)
            absM[i] = std::abs(M[i]);

        
        tol = *max_element(absM.begin(), absM.end());

        //cout << tol << endl;

        iter+=1;
    }    

    // Loadflow: Bus power injections, line & power flows

    unsigned int nl=fb.size();
    
    vc Vm;
    for(i=0;i<V.size();i++)
        Vm.push_back(V[i]*cos(del[i]) + V[i]*sin(del[i])*1i);

    for (auto& f : del) { 
        f = 180/M_PI*f;
    }

    matr Zeta(size_z,vd(10,0));
    vector<vc> Iij(N,vc(N,0)), Sij(N,vc(N,0));;

    // Bus current injections
    vc I;
    vd Im,Ia;

    cvecmultiply(ybus,Vm,I);

    for (i = 0; i < I.size(); i++) {
        Im.push_back(std::abs(I[i]));
        Ia.push_back(std::arg(I[i]));
    }

    // Line Current Flows
    for (i=0; i<nl; i++) {
        Iij[int(fb[i])-1][int(tb[i])-1] = - (Vm[int(fb[i])-1] - Vm[int(tb[i])-1])*ybus[int(fb[i])-1][int(tb[i])-1];
        Iij[int(tb[i])-1][int(fb[i])-1] = - Iij[int(fb[i])-1][int(tb[i])-1];
    }

    // Line Power Flows
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            if (i != j) {
                Sij[i][j] = Vm[i]*(std::conj(Iij[i][j]))*cbMVA; 
            }
        }
    }

    m = 0;
    n = 0;
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            if (j > i) {
                if (Sij[i][j] != (0. + 0i)) {
                    Zeta[m][0] = i+1;
                    Zeta[m][1] = j+1;
                    Zeta[m][2] = Sij[i][j].real();
                    Zeta[m][3] = Sij[i][j].imag();
                    m+=1;
                }
            }    
            else if (i > j) {
                if (Sij[i][j] != (0. + 0i)) {
                    Zeta[n][4] = i+1;
                    Zeta[n][5] = j+1;
                    Zeta[n][6] = Sij[i][j].real();
                    Zeta[n][7] = Sij[i][j].imag();
                    n+=1;
                }
            }    
        }
    }
    
    // Line losses
    vc Lij;
    for (i=0; i<nl; i++) {
        Lij.push_back(Sij[int(fb[i])-1][int(tb[i])-1] + Sij[int(tb[i])-1][int(fb[i])-1]);
        Zeta[i][8] = Lij[i].real();
        Zeta[i][9] = Lij[i].imag();
    }

    // Bus power injections
    vc Si(N,0);
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            Si[i] += (std::conj(Vm[i]))*Vm[j]*ybus[i][j]*cbMVA; 
        }
    }

    vd Pi, Qi;
    for (i=0; i<N; i++) {
        Pi.push_back(Si[i].real());
        Qi.push_back(-Si[i].imag());
        Pl[i]*=bMVA;
        Ql[i]*=bMVA;
    }

    add(Pi, Pl, Pg);
    add(Qi, Ql, Qg);

    double Lijtotalp = 0, Lijtotalq = 0;
    for (i=0; i<Zeta.size(); i++) {
        Lijtotalp += Zeta[i][8];
        Lijtotalq += Zeta[i][9];
    }
    cout << "#######################################################################################################" << endl;
    cout << "-------------------------------------------------------------------------------------------------------" << endl;
    cout << "--------------------------------Newton Raphson Loadflow Analysis---------------------------------------" << endl;
    cout << "-------------------------------------------------------------------------------------------------------" << endl;
    cout << "| Bus |    V   |   Angle    |          Injection         |          Generation       |      Load      |" << endl;
    cout << "| No  |   pu   |   Degree   |      MW     |     MVar     |     MW      |     MVar    |   MW   |  MVar |" << endl;
    cout << "-------------------------------------------------------------------------------------------------------" << endl;
    for(i=0; i<N; i++) {
        cout << std::setw(4) << i+1 << std::setw(11) << V[i] << std::setw(12) << del[i] << std::setw(14) << Pi[i];
        cout << std::setw(15) << Qi[i] << std::setw(14) << Pg[i] << std::setw(14) << Qg[i];
        cout << std::setw(8) << Pl[i] << std::setw(8) << Ql[i] << "\n";
    }
    cout << "-------------------------------------------------------------------------------------------------------" << endl;
    cout << "Total" << std::setw(36) << sumVector(Pi) << std::setw(15) << sumVector(Qi) << std::setw(14);
    cout << sumVector(Pg) << std::setw(14) << sumVector(Qg) << std::setw(8);
    cout << sumVector(Pl) << std::setw(8) << sumVector(Ql) << "\n";
    cout << endl;
    cout << "-------------------------------------------------------------------------------------------------------" << endl;
    cout << "--------------------------------------Line Flow and Losses---------------------------------------------" << endl;
    cout << "-------------------------------------------------------------------------------------------------------" << endl;
    cout << "| From | To  |     P      |     Q     | From | To  |     P      |    Q      |        Line Loss        |" << endl;
    cout << "| Bus  | Bus |    MW      |    MVar   | Bus  | Bus |     MW     |   MVar    |     MW      |    MVar   |" << endl;
    cout << "-------------------------------------------------------------------------------------------------------" << endl;
    for(i=0; i<nl; i++) {
        cout << std::setw(4) << Zeta[i][0] << std::setw(7) << Zeta[i][1] << std::setw(14) << Zeta[i][2] << std::setw(12) << Zeta[i][3];
        cout << std::setw(5) << Zeta[i][4] << std::setw(7) << Zeta[i][5] << std::setw(14) << Zeta[i][6];
        cout << std::setw(12) << Zeta[i][7] << std::setw(14) << Zeta[i][8] << std::setw(12) << Zeta[i][9] << "\n";
    }
    cout << "-------------------------------------------------------------------------------------------------------" << endl;
    cout << "Total" << std::setw(84) << Lijtotalp << std::setw(12) << Lijtotalq << "\n";
    cout << "-------------------------------------------------------------------------------------------------------" << endl;
    cout << "#######################################################################################################" << endl;








        





    return 0;
}
catch (const char * e)
{
    cout << "Error: " << e << '\n';
    return 1;
}
catch (string e)
{
    cout << "Error: " << e << '\n';
    return 1;
}
catch (...)
{
    cout << "Some Error\n";
    return 1;
}

