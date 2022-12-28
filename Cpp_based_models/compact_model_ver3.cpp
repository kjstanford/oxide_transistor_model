// created as a precursor to va implementation; 
// based on ITO_analytical_model_ver2.py (Dual Exponential Trap Distribution & self written NR solver);
// Boole's rule for integration

#include <bits/stdc++.h>
using namespace std;

// set PATH=%PATH%;C:\msys64\mingw64\bin
// g++ -std=c++20 compact_model_ver3.cpp -o compact
// compact

double sgamma(double x)
{
    double _a[30] = { 1.00000000000000000000, 0.57721566490153286061, -0.65587807152025388108,
        -0.04200263503409523553, 0.16653861138229148950, -0.04219773455554433675,
        -0.00962197152787697356, 0.00721894324666309954, -0.00116516759185906511,
        -0.00021524167411495097, 0.00012805028238811619, -0.00002013485478078824,
        -0.00000125049348214267, 0.00000113302723198170, -0.00000020563384169776,
        0.00000000611609510448, 0.00000000500200764447, -0.00000000118127457049,
        0.00000000010434267117, 0.00000000000778226344, -0.00000000000369680562,
        0.00000000000051003703, -0.00000000000002058326, -0.00000000000000534812,
        0.00000000000000122678, -0.00000000000000011813, 0.00000000000000000119,
        0.00000000000000000141, -0.00000000000000000023, 0.00000000000000000002
            };
    double y = x - 1.0;
    double sm = _a[29];
    for (int i=28; i >= 0; i--)
    {
        sm = sm * y + _a[i];
    }
    return 1.0 / sm;
}

// Note convergence restrictions: abs(x) < 1 and c not a negative integer or zero
double hyp21_series( double a, double b, double c, double x )
{
    double TOLERANCE = 1.0e-10;
    double term = a * b * x / c;
    double value = 1.0 + term;
    long long int n = 1;

    while ( abs( term ) > TOLERANCE )
    {
        a++, b++, c++, n++;
        term *= a * b * x / (c * n);
        value += term;
    }

    return value;
}

double hyp21_self( double a, double b, double c, double x )
{
    double value;
    if (x > 0)
        throw std::invalid_argument("can't solve for positive x");
    if (abs(x) < 1)
    {
        // cout << "abs(x) < 1" << endl;
        value = hyp21_series(a,b,c,x);
    }
    if (abs(x) > 1)
    {
        // cout << "abs(x) > 1" << endl;
        value = pow(-x,-a)*(sgamma(c)*sgamma(b-a))/(sgamma(b)*sgamma(c-a))*hyp21_series(a,a-c+1,a-b+1,1/x);
        value += pow(-x,-b)*(sgamma(c)*sgamma(a-b))/(sgamma(a)*sgamma(c-b))*hyp21_series(b-c+1,b,b-a+1,1/x);
    }
    if (abs(x) == 1)
    {
        cout << "abs(x) = 1" << endl;
        // throw std::invalid_argument("can't solve for positive x");
        value = pow(1-x,-a)*hyp21_series(a,c-b,c,x/(x-1));
    }
    return value;
}

double NR_fsolve(auto f, auto df, double z0)
{
    double z = z0;
    double TOLERANCE = 1.0e-10;
    double ALPHA = 0.5;
    double zprev = DBL_MAX;
    int c = 0;
    while ( abs(z-zprev) > TOLERANCE && c < 200)
    {
        zprev = z;
        z = zprev - ALPHA*f(zprev)/df(zprev);
        c++;
    }
    if ( abs(z-zprev) > TOLERANCE )
    {
        throw std::invalid_argument("Could not solve :(");
    }
    return z;
}

vector<double> ITO_compact_model(double Vgs, double Vds) {
    // constants
    double kB = 8.617e-5; //eV/K
    double q = 1.6e-19; //C
    double mo = 9.11e-31; //kg
    double hred = 1.054571817e-34;
    double pi = 3.14159265358979323846;
    double h = hred*2*pi;
    double epso = 8.85e-12;

    // parameters
    double mu_band = 33.5*1e-4;
    double Ntraps = 5.5e16;
    double Ttraps = 400;
    double Ndeep = 5.5e16;
    double Tdeep = 1000;
    double Nch = 1e17;
    double T = 300;
    double tITO = 4.5e-9;
    double tDE = 5.3e-9;
    double me = 0.3*mo;
    double kDE = 16;
    double kITO = 9;
    double phiM = 5.2;
    double chiS = 4.3;
    double L = 1;
    double W = 1;

    double Cox = epso*kDE/tDE;
    double g2D = (me*q)/(pi*pow(hred,2));
    double Etraps = kB*Ttraps;
    double E1 = pow(hred*pi,2)/(2*q*me*pow(tITO,2));
    double E2 = pow(hred*pi*2,2)/(2*q*me*pow(tITO,2));
    double E3 = pow(hred*pi*3,2)/(2*q*me*pow(tITO,2));

    auto gfree = [=] (double E){
        return g2D*((E>=E1)+(E>=E2)+(E>=E3));
    };
    auto gtraps = [=] (double E){
        return (Ntraps/Etraps)*exp((E-E1)/Etraps)*(E<E1);
    };

    auto n2D = [=] (double delE){
        return g2D*kB*T*log(1+exp(-(delE)/(kB*T)));
    };
    auto dn2D = [=] (double delE){
        return g2D*exp(-(delE)/(kB*T))/(1+exp(-(delE)/(kB*T)));
    };
    auto nfree = [=] (double Ef, double phi){
        return n2D(E1-phi-Ef) + n2D(E2-phi-Ef) + n2D(E3-phi-Ef);
    };
    auto dnfree = [=] (double Ef, double phi){
        return dn2D(E1-phi-Ef) + dn2D(E2-phi-Ef) + dn2D(E3-phi-Ef);
    };

    auto ntraps = [=] (double Ef, double phi){
        return Ntraps*hyp21_self(1,T/Ttraps,1+T/Ttraps,-exp((E1-phi-Ef)/(kB*T)));
    };
    auto dntraps = [=] (double Ef, double phi){
        return (Ntraps/(kB*T))*exp((E1-phi-Ef)/(kB*T))*
            hyp21_self(2,1+T/Ttraps,2+T/Ttraps,-exp((E1-phi-Ef)/(kB*T))/(kB*T*(Ttraps/T)*(1+T/Ttraps)));
    };

    auto ndeep = [=] (double Ef, double phi){
        return Ndeep*hyp21_self(1,T/Tdeep,1+T/Tdeep,-exp((E1-phi-Ef)/(kB*T)));
    };
    auto dndeep = [=] (double Ef, double phi){
        return (Ndeep/(kB*T))*exp((E1-phi-Ef)/(kB*T))*
            hyp21_self(2,1+T/Tdeep,2+T/Tdeep,-exp((E1-phi-Ef)/(kB*T))/(kB*T*(Tdeep/T)*(1+T/Tdeep)));
    };
    
    auto f_Ef_eqn = [=] (double Ef){
        return (nfree(Ef,0)+ntraps(Ef,0)+ndeep(Ef,0)-Nch)/Nch;
    };
    auto df_Ef_eqn = [=] (double Ef){
        return (dnfree(Ef,0)+dntraps(Ef,0)+dndeep(Ef,0))/Nch;
    }; 

    double Efermi = NR_fsolve(f_Ef_eqn, df_Ef_eqn, 0);
    cout << Efermi << endl;
    double phiMS = phiM - (chiS-Efermi);
    
    double V = 0;
    auto f = [=] (double phi){
        return Cox*(Vgs-phiMS-phi) - q*(nfree(Efermi,phi-V)+ntraps(Efermi,phi-V)+ndeep(Efermi,phi-V)-Nch);
    };
    auto df = [=] (double phi){
        return -Cox - q*(dnfree(Efermi,phi-V)+dntraps(Efermi,phi-V)+dndeep(Efermi,phi-V));
    };
    double phi_fn = NR_fsolve(f, df, 0);
    auto Id_fn = [=] (double V){
        return (W/L)*q*mu_band*nfree(Efermi,phi_fn-V);
    };
    double IdV_prev = Id_fn(V);

    vector<double> rval (2,0);
    rval[0] = phi_fn;

    double Id_return = 0;
    int N = 100;
    const int NNC = 4;
    for (int i = 0; i < N; i+=NNC) {
        double IdV[NNC+1];
        IdV[0] = IdV_prev;
        for (int j=1; j <= NNC; j++) {
            double V = (i+j)*Vds/N;
            auto f = [=] (double phi){
                return Cox*(Vgs-phiMS-phi) - q*(nfree(Efermi,phi-V)+ntraps(Efermi,phi-V)+ndeep(Efermi,phi-V)-Nch);
            };
            auto df = [=] (double phi){
                return -Cox - q*(dnfree(Efermi,phi-V)+dntraps(Efermi,phi-V)+dndeep(Efermi,phi-V));
            };
            double phi_fn = NR_fsolve(f, df, 0);
            auto Id_fn = [=] (double V){
                return (W/L)*q*mu_band*nfree(Efermi,phi_fn-V);
            };
            IdV[j] = Id_fn(V);  
        }
        IdV_prev = IdV[NNC];
        // cout << IdV_prev << endl;
        double qty = (2.0/45.0)*(Vds/N)*(7.0*IdV[0]+32.0*IdV[1]+12.0*IdV[2]+32.0*IdV[3]+7.0*IdV[4]);
        Id_return += qty; 
    }

    rval[1] = Id_return;
    return rval;
}

int main() {
    double Vgs, Vds;
    cin >> Vgs;
    cin >> Vds;
    vector<double> a = ITO_compact_model(Vgs,Vds);
    cout.precision(8);
    cout << a[0] << " " << a[1];
    return 0;
}