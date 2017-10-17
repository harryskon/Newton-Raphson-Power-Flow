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
void zdatas(int num, vd zd[], vd out[], matr G, matr B, matr bbus);

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

void printMatrix(matr A) 
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

    if (N==14)
        sn = &s14;
    else if (N==30)
        sn = &s30;
    else
        throw "Bad input IEEE system argument";

    
    //Form Admittance (Y) And Impedance (Z) (inv(ybus)) Bus Formation
    unsigned int i,j;
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
    unsigned int QG;

    while (iter < 5) {
        
        vd P(N,0), Q(N,0);

        // Calculate P and Q
        for(i=0;i<N;i++) {
            for(j=0;j<N;j++){
                P[i] = P[i] + V[i]*V[j]*(G[i][j]*cos(del[i]-del[j]) + B[i][j]*sin(del[i]-del[j]));
                Q[i] = Q[i] + V[i]*V[j]*(G[i][j]*sin(del[i]-del[j]) - B[i][j]*cos(del[i]-del[j]));
            }
        }

        // Check Q-limit violations
        if ((iter > 2) && (iter < 8)) {
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



        iter+=1;
    }    



    //vd zd[6], out[2];


    // for(unsigned int runs=0; runs<nruns; runs++) {

    //     vd V(N,1), del(N,0);
    //     vd E,vi,pi,qi,pf,qf,typez,z,fbus,tbus,rii,idxm,hx;
        
    //     E.insert(E.end(), del.begin()+1, del.end());
    //     E.insert(E.end(), V.begin(), V.end());
        
    //     zdatas(N,zd,out,G,B,bbus);
    //     for(i=0; i<2; i++){
    //         out[i].clear();
    //     }
        
    //     /* Traditional Measurement Data..
    //       Type: Vi - 1, Pi - 2, Qi - 3, Pij - 4, Qij - 5;
    //      |Msnt |Type | Value | From | To | Rii | */
    
    //     for(i=0; i<zd[0].size(); i++ ) {
    //         idxm.push_back(zd[0][i]);
    //         typez.push_back(zd[1][i]);
    //         z.push_back(zd[2][i]);
    //         fbus.push_back(zd[3][i]);
    //         tbus.push_back(zd[4][i]);
    //         rii.push_back(zd[5][i]);
    //         switch (int(zd[1][i])) {
    //             case 1: vi.push_back(zd[0][i]); break;
    //             case 2: pi.push_back(zd[0][i]); break;
    //             case 3: qi.push_back(zd[0][i]); break;
    //             case 4: pf.push_back(zd[0][i]); break;
    //             case 5: qf.push_back(zd[0][i]); break;
    //             default: throw "measurement input file is corrupted";

    //         }
    //     }

    //     matr diagrii(rii.size(),vd(rii.size()));
    //     matr invdiagrii(rii.size(),vd(rii.size()));

    //     //diagonal matrix of rii
    //     for(i=0; i<rii.size(); i++ ){
    //         for(j=0; j<rii.size(); j++ ) {
    //             if (i == j) {
    //                 diagrii[i][j] = rii[i];
    //                 invdiagrii[i][j] = 1/rii[i];
    //             }
    //             else {
    //                 diagrii[i][j] = 0;
    //                 invdiagrii[i][j] = 0;
    //             }
    //         }    
    //     }

    //     int iter = 1;
    //     int x;
    //     double tol = 5.0;
    //     vd h[5];
    //     vd residue(z.size());
    //     double J = 0;

    //     while (tol > 1e-4) {

    //         //Measurement function h
    //         for(i=0; i<vi.size(); i++){
    //             x = vi[i]-1;
    //             if (x < 0) x = 0;
    //             m = fbus[x] - 1;
    //             h[0].push_back(V[m]);
    //         }
    //         h[1] = vd(pi.size(),0);
    //         h[2] = vd(qi.size(),0);
    //         h[3] = vd(pf.size(),0);
    //         h[4] = vd(qf.size(),0);
            
    //         for(i=0; i<pi.size(); i++) {
    //             x = pi[i]-1;
    //             if (x < 0) x = 0;
    //             m = fbus[x] - 1;
    //             for(j=0; j<N; j++) {
    //                 h[1][i] = h[1][i] + V[m]*V[j]*(G[m][j]*cos(del[m] - del[j]) + B[m][j]*sin(del[m] - del[j]));
    //             }
    //         }

    //         for(i=0; i<qi.size(); i++) {
    //             x = qi[i]-1;
    //             if (x < 0) x = 0;
    //             m = fbus[x] - 1;
    //             if (m < 0) m = 0;
    //             for(j=0; j<N; j++) {
    //                 h[2][i] = h[2][i] + V[m]*V[j]*(G[m][j]*sin(del[m] - del[j]) - B[m][j]*cos(del[m] - del[j]));
    //             }
    //         }

    //         for(i=0; i<pf.size(); i++) {
    //             x = pf[i]-1;
    //             if (x < 0) x = 0;
    //             m = fbus[x] - 1;
    //             n = tbus[x] - 1;
    //             if (n < 0) n = 0;
    //             if (m < 0) m = 0;
    //             h[3][i] = -V[m]*V[m]*G[m][n] - V[m]*V[n]*(-G[m][n]*cos(del[m] - del[n]) - B[m][n]*sin(del[m] - del[n]));
    //         }
  
    //         for(i=0; i<qf.size(); i++) {
    //             x = qf[i]-1;
    //             if (x < 0) x = 0;
    //             m = fbus[x] - 1;
    //             n = tbus[x] - 1;
    //             if (n < 0) n = 0;
    //             if (m < 0) m = 0;
    //             h[4][i] = -(V[m]*V[m])*(-B[m][n] + bbus[m][n]) - V[m]*V[n]*(-G[m][n]*sin(del[m] - del[n]) + B[m][n]*cos(del[m] - del[n]));
    //         }
            
    //         for(i=0; i<5; i++){
    //             hx.insert(hx.end(), h[i].begin(), h[i].end());
    //         }
            
    //         //get residue
    //         subtract(z,hx,residue);

    //         // Jacobian..
    //         // H11 - Derivative of V with respect to angles.. All Zeros
    //         matr H11(vi.size(),vd(N-1,0));
    //         matr H12(vi.size(),vd(N,0));
    //         matr H21(pi.size(),vd(N-1,0));
    //         matr H22(pi.size(),vd(N,0));
    //         matr H31(qi.size(),vd(N-1,0));
    //         matr H32(qi.size(),vd(N,0));
    //         matr H41(pf.size(),vd(N-1,0));
    //         matr H42(pf.size(),vd(N,0));
    //         matr H51(qf.size(),vd(N-1,0));
    //         matr H52(qf.size(),vd(N,0));

    //         i = H11.size()+H21.size()+H31.size()+H41.size()+H51.size();            
    //         matr H(i);
            
    //         //H12 - Derivative of V with respect to V..
    //         for(i=0;i<vi.size();i++) {
    //             for(j=0;j<N;j++){
    //                 if (i==j) {
    //                     H12[i][j] = 1;
    //                 }    
    //             }
    //         }

    //        // H21 - Derivative of Real Power Injections with Angles..
    //         for(i=0;i<pi.size();i++) {
    //             x = pi[i]-1;
    //             if (x < 0) x = 0;
    //             m = fbus[x] - 1;
    //             for(j=0;j<N-1;j++){
    //                 if (j+1 == m) {
    //                     for(k=0;k<N;k++){
    //                         H21[i][j] = H21[i][j] + V[m]*V[k]*(-G[m][k]*sin(del[m]-del[k]) + B[m][k]*cos(del[m]-del[k]));
    //                     }
    //                     H21[i][j] = H21[i][j] - V[m]*V[m]*B[m][m];  
    //                 }
    //                 else {
    //                     H21[i][j] = V[m]*V[j+1]*(G[m][j+1]*sin(del[m]-del[j+1]) - B[m][j+1]*cos(del[m]-del[j+1]));
    //                 }    
    //             }
    //         }


    //         // H22 - Derivative of Real Power Injections with V..
    //         for(i=0;i<pi.size();i++) {
    //             x = pi[i]-1;
    //             if (x < 0) x = 0;
    //             m = fbus[x] - 1;
    //             for(j=0;j<N;j++){
    //                 if (j==m) {
    //                     for(k=0;k<N;k++){
    //                         H22[i][j] = H22[i][j] + V[k]*(G[m][k]*cos(del[m]-del[k]) + B[m][k]*sin(del[m]-del[k]));
    //                     }
    //                     H22[i][j] = H22[i][j] + V[m]*G[m][m];  
    //                 }
    //                 else {
    //                     H22[i][j] = V[m]*(G[m][j]*cos(del[m]-del[j]) + B[m][j]*sin(del[m]-del[j]));
    //                 }    
    //             }
    //         }

            
    //         // H31 - Derivative of Reactive Power Injections with Angles..
    //         for(i=0;i<qi.size();i++) {
    //             x = qi[i]-1;
    //             if (x < 0) x = 0;
    //             m = fbus[x] - 1;
    //             for(j=0;j<N-1;j++){
    //                 if (j+1 == m) {
    //                     for(k=0;k<N;k++){
    //                         H31[i][j] = H31[i][j] + V[m]*V[k]*(G[m][k]*cos(del[m]-del[k]) + B[m][k]*sin(del[m]-del[k]));
    //                     }
    //                     H31[i][j] = H31[i][j] - V[m]*V[m]*G[m][m];  
    //                 }
    //                 else {
    //                     H31[i][j] = V[m]*V[j+1]*(-G[m][j+1]*cos(del[m]-del[j+1]) - B[m][j+1]*sin(del[m]-del[j+1]));
    //                 }    
    //             }
    //         }
            
    //         // H32 - Derivative of Reactive Power Injections with V..
    //         for(i=0;i<qi.size();i++) {
    //             x = qi[i]-1;
    //             if (x < 0) x = 0;
    //             m = fbus[x] - 1;
    //             for(j=0;j<N;j++){
    //                 if (j==m) {
    //                     for(k=0;k<N;k++){
    //                         H32[i][j] = H32[i][j] + V[k]*(G[m][k]*sin(del[m]-del[k]) - B[m][k]*cos(del[m]-del[k]));
    //                     }
    //                     H32[i][j] = H32[i][j] - V[m]*B[m][m];  
    //                 }
    //                 else {
    //                     H32[i][j] = V[m]*(G[m][j]*sin(del[m]-del[j]) - B[m][j]*cos(del[m]-del[j]));
    //                 }    
    //             }
    //         }

    //         // H41 - Derivative of Real Power Flows with Angles..
    //         for(i=0;i<pf.size();i++) {
    //             x = pf[i]-1;
    //             if (x < 0) x = 0;
    //             m = fbus[x] - 1;
    //             n = tbus[x] - 1;
    //             if (m < 0) m = 0;
    //             if (n < 0) n = 0;
    //             for(j=0;j<N-1;j++){
    //                 if (j+1 == m) {
    //                     H41[i][j] = V[m]*V[n]*(-G[m][n]*sin(del[m]-del[n]) + B[m][n]*cos(del[m]-del[n]));
    //                 }
    //                 else if (j+1 == n){
    //                     H41[i][j] = -V[m]*V[n]*(-G[m][n]*sin(del[m]-del[n]) + B[m][n]*cos(del[m]-del[n]));
    //                 }
    //                 else {
    //                     H41[i][j] = 0;
    //                 }    
    //             }
    //         }

            
    //         // H42 - Derivative of Real Power Flows with V..
    //         for(i=0;i<pf.size();i++) {
    //             x = pf[i]-1;
    //             if (x < 0) x = 0;
    //             m = fbus[x] - 1;
    //             n = tbus[x] - 1;
    //             if (m < 0) m = 0;
    //             if (n < 0) n = 0;
    //             for(j=0;j<N;j++){
    //                 if (j == m) {
    //                     H42[i][j] = -V[n]*(-G[m][n]*cos(del[m]-del[n]) - B[m][n]*sin(del[m]-del[n])) - 2*G[m][n]*V[m];
    //                 }
    //                 else if (j == n){
    //                     H42[i][j] = -V[m]*(-G[m][n]*cos(del[m]-del[n]) - B[m][n]*sin(del[m]-del[n]));
    //                 }
    //                 else {
    //                     H42[i][j] = 0;
    //                 }    
    //             }
    //         }

    //         // H51 - Derivative of Reactive Power Flows with Angles..
    //         for(i=0;i<qf.size();i++) {
    //             x = qf[i]-1;
    //             if (x < 0) x = 0;
    //             m = fbus[x] - 1;
    //             n = tbus[x] - 1;
    //             if (m < 0) m = 0;
    //             if (n < 0) n = 0;
    //             for(j=0;j<N-1;j++){
    //                 if (j+1 == m) {
    //                     H51[i][j] = -V[m]*V[n]*(-G[m][n]*cos(del[m]-del[n]) - B[m][n]*sin(del[m]-del[n]));
    //                 }
    //                 else if (j+1 == n){
    //                     H51[i][j] = V[m]*V[n]*(-G[m][n]*cos(del[m]-del[n]) - B[m][n]*sin(del[m]-del[n]));
    //                 }
    //                 else {
    //                     H51[i][j] = 0;
    //                 }    
    //             }
    //         }
            
    //         // H52 - Derivative of Reactive Power Flows with V..
    //         for(i=0;i<qf.size();i++) {
    //             x = qf[i]-1;
    //             if (x < 0) x = 0;
    //             m = fbus[x] - 1;
    //             n = tbus[x] - 1;
    //             if (m < 0) m = 0;
    //             if (n < 0) n = 0;
    //             for(j=0;j<N;j++){
    //                 if (j == m) {
    //                     H52[i][j] = -V[n]*(-G[m][n]*sin(del[m]-del[n]) + B[m][n]*cos(del[m]-del[n])) - 2*V[m]*(-B[m][n] + bbus[m][n]);
    //                 }
    //                 else if (j == n){
    //                     H52[i][j] = -V[m]*(-G[m][n]*sin(del[m]-del[n]) + B[m][n]*cos(del[m]-del[n]));
    //                 }
    //                 else {
    //                     H52[i][j] = 0;
    //                 }    
    //             }
    //         }

    //         //Measurement Jacobian H
    //         for(i=0;i<H11.size();i++) {
    //             H[i].insert(H[i].end(), H11[i].begin(), H11[i].end());
    //             H[i].insert(H[i].end(), H12[i].begin(), H12[i].end());
    //         }
    //         int tmp = H11.size();
    //         for(i=0;i<H21.size();i++) {
    //             H[tmp+i].insert(H[tmp+i].end(), H21[i].begin(), H21[i].end());
    //             H[tmp+i].insert(H[tmp+i].end(), H22[i].begin(), H22[i].end());
    //         }
    //         tmp = H11.size() + H21.size();
    //         for(i=0;i<H31.size();i++) {
    //             H[tmp+i].insert(H[tmp+i].end(), H31[i].begin(), H31[i].end());
    //             H[tmp+i].insert(H[tmp+i].end(), H32[i].begin(), H32[i].end());
    //         }
    //         tmp = H11.size() + H21.size() + H31.size();
    //         for(i=0;i<H41.size();i++) {
    //             H[tmp+i].insert(H[tmp+i].end(), H41[i].begin(), H41[i].end());
    //             H[tmp+i].insert(H[tmp+i].end(), H42[i].begin(), H42[i].end());
    //         }
    //         tmp = H11.size() + H21.size() + H31.size() + H41.size();
    //         for(i=0;i<H51.size();i++) {
    //             H[tmp+i].insert(H[tmp+i].end(), H51[i].begin(), H51[i].end());
    //             H[tmp+i].insert(H[tmp+i].end(), H52[i].begin(), H52[i].end());
    //         }
            

    //          //Gm = H'*inv(diagrii)*H; Gain Matrix, Gm..
    //         i=H[0].size();
    //         j=H.size();
    //         matr Htrans(i,vd(j));
            
    //         i=Htrans.size();
    //         j=invdiagrii[0].size();
    //         matr temp(i,vd(j));

    //         i=Htrans.size();
    //         j=H[0].size();
    //         matr Gm(i,vd(j));

    //         for (i = 0; i < H[0].size(); i++) {
    //             for (j = 0; j < H.size(); j++) {
    //                 Htrans[i][j] = H[j][i];
    //             }
    //         }

    //         for(i = 0; i < temp.size(); i++)
    //             for(j = 0; j < temp[0].size(); j++)
    //                temp[i][j] =  matrmultiply(Htrans, invdiagrii, i, j);

    //         for(i = 0; i < Gm.size(); i++)
    //             for(j = 0; j < Gm[0].size(); j++)
    //                Gm[i][j] =  matrmultiply(temp, H, i, j);

    //         unsigned int size_Gm = Gm.size();

    //         // Covariance matrix (not used now)
    //         vd CvE(size_Gm);
    //         for(i=0; i<Gm.size(); i++ ){
    //             for(j=0; j<Gm[0].size(); j++ ) {
    //                 if (i == j) {
    //                     CvE[i] = Gm[i][j];
    //                 }
    //             }    
    //         }

    //         /*Objective Function: J = sum(inv(Ri)*r.^2) (not used now) */
    //         vd rsq(residue.size());
    //         vd t;
    //         multiply(residue,residue,rsq);

    //         vecmultiply(invdiagrii,rsq,t);
    //         std::for_each(t.begin(), t.end(), [&] (double n) { J += n;});
        
    //         // State vector
    //         t.clear();
    //         vecmultiply(temp,residue,t);

    //         // calc inv(Gm)
    //         Matrix RealA2(size_Gm,size_Gm);
    //         Matrix ImagA2(size_Gm,size_Gm);
    //         Matrix RealA2inv(size_Gm,size_Gm);
    //         Matrix ImagA2inv(size_Gm,size_Gm);

    //         for (i=0; i<size_Gm; i++){
    //             for (j=0; j<size_Gm; j++){
    //                 RealA2(i+1,j+1) = Gm[i][j];
    //                 ImagA2(i+1,j+1) = 0;
    //             }   
    //         }
    //         cinv(RealA2, ImagA2, RealA2inv, ImagA2inv);

    //         matr invGm(size_Gm,vd(size_Gm));

    //         //inv(Gm)
    //         for (i=0; i<size_Gm; i++){
    //             for (j=0; j<size_Gm; j++){
    //                 invGm[i][j] = RealA2inv(i+1,j+1);
    //             }   
    //         }

    //         //dE = inv(Gm)*(H'*inv(Ri)*r)
    //         vd dE;
    //         vecmultiply(invGm,t,dE);

    //         // E = E + dE
    //         add(E,dE,E);

    //         //del(2:end) = E(1:nbus-1)
    //         //V = E(nbus:end)   
    //         del.erase(del.begin()+1, del.end());
    //         V.clear();
    //         del.insert(del.end(), E.begin(), E.begin()+N-1);
    //         V.insert(V.end(), E.begin()+N-1, E.end());

    //         iter = iter + 1;

    //         //tol = max(abs(dE))
    //         vd absdE(dE.size());
    //         for (i = 0; i < dE.size(); i++){
    //             if (dE[i] < 0) 
    //                 absdE[i] = (-1)*dE[i];
    //             else
    //                 absdE[i] = dE[i];
    //         }
    //         tol = *max_element(absdE.begin(), absdE.end());  

    //         for(i=0; i<5; i++)
    //             h[i].clear();

    //         hx.clear();          
    //     }
        
    //     out[0].insert(out[0].end(), V.begin(), V.end());
    //     out[1].insert(out[1].end(), del.begin(), del.end());

    //     for (auto& f : del) { f = 180/M_PI*f;}

    //     cout << "-------- State Estimation ------------------" << endl;
    //     cout << "--------------------------" << endl;
    //     cout << "| Bus |    V   |  Angle  | " << endl;
    //     cout << "| No  |   pu   |  Degree | " << endl;
    //     cout << "--------------------------" << endl;
    //     for(i=0; i<N; i++)
    //         cout << string(3,' ') << i+1 << string(3,' ') << V[i] << string(3,' ') << del[i] << "\n";
    //     cout << "---------------------------------------------" << endl;

    //     V.clear();
    //     del.clear();
    //     E.clear();
    // }


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

