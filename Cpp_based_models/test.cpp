#include <bits/stdc++.h>
using namespace std;

// set PATH=%PATH%;C:\msys64\mingw64\bin
// g++ -std=c++20 test.cpp -o test
// test

// vector<double> CalculateCoefficients(int numCoeff)
// {
//     vector<double> c(numCoeff);
//     double k1_factrl = 1.0;
//     c[0] = sqrt(2.0 * M_PI);
//     for(size_t k=1; k < numCoeff; k++)
//     {
//         c[k] = exp(numCoeff-k) * pow(numCoeff-k, k-0.5) / k1_factrl;
//         k1_factrl *= -(double)k;
//     }
//     return c;
// }

// // The Spouge approximation
// double Gamma(const std::vector<double>& coeffs, double x)
// {
//         const size_t numCoeff = coeffs.size();
//         double accm = coeffs[0];
//         for(size_t k=1; k < numCoeff; k++)
//         {
//             accm += coeffs[k] / ( x + k );
//         }
//         accm *= exp(-(x+numCoeff)) * pow(x+numCoeff, x+0.5);
//         return accm/x;
// }

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
    const double TOLERANCE = 1.0e-10;
    double term = a * b * x / c;
    double value = 1.0 + term;
    long long int n = 1;

    while ( abs( term ) > TOLERANCE )
    {
        a++, b++, c++, n++;
        term *= a * b * x / c / n;
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
        // cout << "abs(x) = 1" << endl;
        // throw std::invalid_argument("can't solve for positive x");
        value = pow(1-x,-a)*hyp21_series(a,c-b,c,x/(x-1));
    }
    return value;
}

double NR_fsolve(auto f, auto df, double z0)
{
    double z = z0;
    const double TOLERANCE = 1.0e-10;
    const double ALPHA = 0.8;
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

int main()
{
    double a = 1, b = 0.05;
    double c = a + b;
    double x;
    cin >> x;
    cout << "F( " << a << ", " << b << "; " << c << ", " << x << " ) = ";
    cout << hyp21_self( a, b, c, x ) << '\n';

    auto f = [=] (double z){
        return pow(z,3)-2*z-5;
    };
    auto df = [=] (double z){
        return 3*pow(z,2)-2;
    };
    cout << NR_fsolve(f,df,2);
    return 0;
}






// #include <Python.h>
// #include <stdlib.h>
// int main()
// {
//    // Set PYTHONPATH TO working directory
//    setenv("PYTHONPATH",".",1);

//    PyObject *pName, *pModule, *pDict, *pFunc, *pValue, *presult;


//    // Initialize the Python Interpreter
//    Py_Initialize();


//    // Build the name object
//    pName = PyString_FromString((char*)"arbName");

//    // Load the module object
//    pModule = PyImport_Import(pName);


//    // pDict is a borrowed reference 
//    pDict = PyModule_GetDict(pModule);


//    // pFunc is also a borrowed reference 
//    pFunc = PyDict_GetItemString(pDict, (char*)"someFunction");

//    if (PyCallable_Check(pFunc))
//    {
//        pValue=Py_BuildValue("(z)",(char*)"something");
//        PyErr_Print();
//        printf("Let's give this a shot!\n");
//        presult=PyObject_CallObject(pFunc,pValue);
//        PyErr_Print();
//    } else 
//    {
//        PyErr_Print();
//    }
//    printf("Result is %d\n",PyInt_AsLong(presult));
//    Py_DECREF(pValue);

//    // Clean up
//    Py_DECREF(pModule);
//    Py_DECREF(pName);

//    // Finish the Python Interpreter
//    Py_Finalize();


//     return 0;
// }