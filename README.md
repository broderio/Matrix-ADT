# Matrix-ADT

Matrix is a C++ library for linear algebra applications.

## Usage

```C++
#include "Matrix.h"

matrix mat_1();     // initializes an empty matrix
matrix mat_2(3);    // initializes an empty 3x3 matrix
matrix mat_3(3,2);  // initializes an empty 3x2 matrix

double a[9] = {1, 2, 3,
               4, 5, 6,
               7, 8, 9};
               
double b[6] = {1, 2,
               3, 4,
               5, 6};

matrix mat_4(3,a);    // initializes a 3x3 matrix with data
matrix mat_5(3,2,b);  // initializes a 3x2 matrix with data

matrix L(3), U(3);
mat_4.LU(L,U);      // Factors mat_4 into lower and upper triangular matrices

matrix Q(3), R(3);
mat_4.QR(Q,R);      // Factors mat_4 into an orthogonal matrix Q and upper triangular matrix R

double A_data[9] = {31, 25, 13, 64, 53, 46, 97, 68, 19};
double b_data[3] = {1, 2, 3};
matrix A(3,A_data), b(3,1,b_data), x;
x = linearSolve(A,b);   // Returns the solution x to the linear system A*x = b

double x1[5] = {1, 2, 4, 5, 7};
double y1[5] = {4, 8, 10, 12, 18};
matrix X1(5,1,x1), Y1(5,1,y1);
matrix C = fit(X1,Y1);  // Returns a 2x1 matrix with coefficients for a linear fit (C1*x + C2)

double x2[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
double y2[11] = {20, 15, 16, 18, 23, 30, 40, 55, 70, 86, 110};
matrix X2(11,1,x2), Y2(11,1,y2);
matrix C = fit(X2,Y2,2);  // Returns a 3x1 matrix with coefficients for a quadratic fit (C1*x^2 + C2*x + C3)
```
