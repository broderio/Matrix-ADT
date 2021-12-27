//
//  Matrix.h
//  Matrix
//
//  Created by Broderick Riopelle on 12/22/21.
//

#ifndef Matrix_h
#define Matrix_h

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include "MatrixClassDef.h"
#include "MatrixErrors.h"

// REQUIRES: N/A
// MODIFIES: N/A
// EFFECTS: Returns true if a is within tolerance of b.
bool is_approx(const double &a, const double &b, const double &tol)
{
    if (a-b > tol) return false;
    return true;
}

// REQUIRES: dim > 0
// MODIFIES: N/A
// EFFECTS:  Returns a dim by dim matrix full of zeros
matrix zeros(const int &dim)
{
    matrix Z(dim);
    for (int i = 0; i < dim; ++i)
    {
        for (int j = 0; j < dim; ++j)
        {
            Z(i,j) = 0;
        }
    }
    return Z;
}

// REQUIRES: row > 0, col > 0
// MODIFIES: N/A
// EFFECTS:  Returns a row by col matrix full of zeros
matrix zeros(const int &row, const int &col)
{
    matrix Z(row,col);
    for (int i = 0; i < row; ++i)
    {
        for (int j = 0; j < col; ++j)
        {
            Z(i,j) = 0;
        }
    }
    return Z;
}

// REQUIRES: dim > 0
// MODIFIES: N/A
// EFFECTS:  Returns a dim by dim matrix full of ones
matrix ones(const int &dim)
{
    matrix Z(dim);
    for (int i = 0; i < dim; ++i)
    {
        for (int j = 0; j < dim; ++j)
        {
            Z(i,j) = 1;
        }
    }
    return Z;
}

// REQUIRES: row > 0, col > 0
// MODIFIES: N/A
// EFFECTS:  Returns a row by col matrix full of ones
matrix ones(const int &row, const int &col)
{
    matrix Z(row,col);
    for (int i = 0; i < row; ++i)
    {
        for (int j = 0; j < col; ++j)
        {
            Z(i,j) = 1;
        }
    }
    return Z;
}

// Constructors

// REQUIRES: N/A
// MODIFIES: Initializes 0 by 0 empty matrix
// EFFECTS: N/A
matrix::matrix()
: _row(0), _col(0), _size(0), data(nullptr), _square(false), _empty(true) {}

// REQUIRES: N/A
// MODIFIES: Initializes empty square matrix with dimension dim
// EFFECTS: N/A
matrix::matrix(const int &dim)
: _row(dim), _col(dim), _size(dim*dim), data(new double[dim*dim]), _square(true),
_empty(false) {}

// REQUIRES: arr must have dim^2 elements
// MODIFIES: Initializes square matrix with dimension dim and with data in arr.
//           The data is assumed to be in order of a 2D matrix starting left to
//           right and going downward.
//           Example: matrix(2, {1, 2, 3, 4}) would result in [ 1 2 ]
//                                                            [ 3 4 ]
// EFFECTS: N/A
matrix::matrix(const int &dim, const double *arr)
: _row(dim), _col(dim), _size(dim*dim), data(new double[dim*dim]), _square(true),
_empty(false)
{
    for (int i = 0; i < size(); ++i) data[i] = arr[i];
}

// REQUIRES: N/A
// MODIFIES: Initializes row_size by col_size empty matrix
// EFFECTS: N/A
matrix::matrix(const int &row_size, const int &col_size)
: _row(row_size), _col(col_size), _size(row_size*col_size),
data(new double[row_size*col_size]), _square(_row == _col), _empty(_size == 0)
{}

// REQUIRES: arr must be of size row_size * col_size
// MODIFIES: Initializes row_size by col_size matrix with data in arr. The data
//           is assumed to be in order of a 2D matrix starting left to right and
//           going downward.
//           Example: matrix(2, 4, {1, 2, 3, 4, 5, 6, 7, 8}) would result in
//                    [ 1 2 3 4 ]
//                    [ 5 6 7 8 ]
// EFFECTS: N/A
matrix::matrix(const int &row_size, const int &col_size, const double *arr)
: _row(row_size), _col(col_size), _size(row_size*col_size),
  data(new double[row_size*col_size]), _square(_row == _col), _empty(_size == 0)
{
      for (int i = 0; i < size(); ++i) data[i] = arr[i];
}


// Big Three

// REQUIRES: N/A
// MODIFIES: Creates a deep copy of other.
// EFFECTS:  N/A
matrix::matrix(const matrix &other)
: _row(other.row()), _col(other.column()), _size(other.size()),
_empty(other.empty()), data(new double[other.size()]), _square(other.is_square())
{
    for (int i = 0; i < row(); ++i)
    {
        for (int j = 0; j < column(); ++j)
        {
            (*this)(i,j) = other.at(i,j);
        }
    }
}

// REQUIRES: N/A
// MODIFIES: Deletes matrix.
// EFFECTS:  N/A
matrix::~matrix()
{
    clear();
}

// REQUIRES: N/A
// MODIFIES: Deletes original matrix and creates a deep copy of rhs.
// EFFECTS:  Returns the matrix after operation.
const matrix &matrix::operator=(const matrix &rhs)
{
    if (this->data == rhs.data) return *this;
    clear();
    _row = rhs.row();
    _col = rhs.column();
    _size = rhs.size();
    _square = rhs.is_square();
    data = new double[_size];
    for (int i = 0; i < row(); ++i)
    {
        for (int j = 0; j < column(); ++j)
        {
            (*this)(i,j) = rhs.at(i,j);
        }
    }
    _empty = rhs.size() > 0;
    return *this;
}


// Overloaded Operators

// REQUIRES: row and col are less than the number of rows and columns in the
//           matrix. row and col are 0 indexed.
// MODIFIES: N/A
// EFFECTS:  Returns a reference to the object in the matrix.
double &matrix::operator()(const int &row, const int &col)
{
    if (row < 0 || row > _row-1 || col < 0 || col > _col-1)
    {
        throw invalid_index("Unable to access elements. Indices out of bounds."
                          " row = "+std::to_string(row)+
                          " col = "+std::to_string(col));
    }
    return data[row*column() + col];
}

// REQUIRES: rhs is the same dimensions as lhs
// MODIFIES: N/A
// EFFECTS:  Returns the sum of two matrices
matrix matrix::operator+(const matrix &rhs) const
{
    if (rhs.column() != column() || rhs.row() != row())
    {
        throw dimension_mismatch("Dimension mismatch. matrix addition"
                              "requires matrices of equal dimensions");
    }
    matrix m(*this);
    for (int i = 0; i < row(); ++i)
    {
        for (int j = 0; j < column(); ++j)
        {
            m(i,j) += rhs.at(i,j);
        }
    }
    return m;
}

// REQUIRES: N/A
// MODIFIES: N/A
// EFFECTS:  Returns the matrix summed with a scalar;
matrix matrix::operator+(const double &a) const
{
    matrix m(*this);
    for (int i = 0; i < row(); ++i)
    {
        for (int j = 0; j < column(); ++j)
        {
            m(i,j) += a;
        }
    }
    return m;
}

// REQUIRES: rhs is the same size and dimension as lhs
// MODIFIES: N/A
// EFFECTS:  Returns the difference of two matrices
matrix matrix::operator-(const matrix &rhs) const
{
    if (rhs.column() != column() || rhs.row() != row())
    {
        throw dimension_mismatch("Dimensions mismatch. matrix subtraction"
                              "requires matrices of equal dimensions");
    }
    matrix m(*this);
    for (int i = 0; i < row(); ++i)
    {
        for (int j = 0; j < column(); ++j)
        {
            m(i,j) -= rhs.at(i,j);
        }
    }
    return m;
}

// REQUIRES: The number of columns in *this is equal to the number of rows in
//           rhs
// MODIFIES: N/A
// EFFECTS:  Returns the dot product of two matrices
matrix matrix::operator*(const matrix &rhs) const
{
    if (rhs.row() != column())
    {
        throw dimension_mismatch("Dimensions mismatch. matrix multiplication "
                              "requires number of columns in first matrix "
                              "to equal number of rows in second matrix.");
    }
    matrix m(_row, rhs.column());
    for (int i = 0; i < _row; ++i)
    {
        for (int j = 0; j < rhs.column(); ++j)
        {
            double sum = 0;
            for (int k = 0; k < _col; ++k)
            {
                sum += this->at(i,k) * rhs.at(k,j);
            }
            m(i,j) = sum;
        }
    }
    return m;
}

// REQUIRES: N/A
// MODIFIES: N/A
// EFFECTS:  Returns the matrix multiplied by a scalar;
matrix matrix::operator*(const double &a) const
{
    matrix m(*this);
    for (int i = 0; i < row(); ++i)
    {
        for (int j = 0; j < column(); ++j)
        {
            m(i,j) *= a;
        }
    }
    return m;
}

// REQUIRES: N/A
// MODIFIES: N/A
// EFFECTS:  Returns the matrix divided by a scalar;
matrix matrix::operator/(const double &a) const
{
    matrix m(*this);
    for (int i = 0; i < row(); ++i)
    {
        for (int j = 0; j < column(); ++j)
        {
            m(i,j) /= a;
        }
    }
    return m;
}

// Functions

// REQUIRES: N/A
// MODIFIES: N/A
// EFFECTS:  Returns number of elements in matrix
int matrix::size() const
{
    return _size;
}

// REQUIRES: N/A
// MODIFIES: N/A
// EFFECTS:  Returns number of rows in matrix
int matrix::row() const
{
    return _row;
}

// REQUIRES: N/A
// MODIFIES: N/A
// EFFECTS:  Returns number of columns in matrix
int matrix::column() const
{
    return _col;
}

// REQUIRES: N/A
// MODIFIES: N/A
// EFFECTS:  Returns true if matrix is square, false otherwise
bool matrix::is_square() const
{
    return _square;
}

// REQUIRES: N/A
// MODIFIES: N/A
// EFFECTS:  Returns true if matrix is empty, false otherwise
bool matrix::empty() const
{
    return _empty;
}

// REQUIRES: N/A
// MODIFIES: N/A
// EFFECTS:  Returns dimension if matrix is square, 0 otherwise
int matrix::dim() const
{
    if (_square) return _row;
    else return 0;
}

// REQUIRES: N/A
// MODIFIES: N/A
// EFFECTS:  Returns true if matrix contains n.
bool matrix::contains(const double &n) const
{
    for (int i = 0; i < row(); ++i)
    {
        for (int j = 0; j < column(); ++j)
        {
            if (is_approx(at(i,j),n)) return true;
        }
    }
    return false;
}

// REQUIRES: matrix must be square. row must be greater than 0 and less than
//           the number of rows and col must be greater than 0 and less than
//           the number of columns.
// MODIFIES: N/A
// EFFECTS:  Returns cofactor matrix of given index.
matrix matrix::cof(const int &row, const int &col) const
{
    if (!is_square())
    {
        throw rectangular_matrix("matrix is not square. Unable to compute "
                                 "cofactor.");
    }
    matrix c(this->row()-1);
    for (int i = 0; i < c.row(); ++i)
    {
        for (int j = 0; j < c.column(); ++j)
        {
            if (i != row && j != col) c(i,j) = at(i,j);
        }
    }
    return c;
}

// REQUIRES: matrix must be square
// MODIFIES: N/A
// EFFECTS:  Returns determinant of matrix
double matrix::det() const
{
    if (!is_square())
    {
        throw rectangular_matrix("matrix is not square. Unable to compute"
                                 "determinant.");
    }
    matrix L(dim()),U(dim());
    double D = 1;
    LU(L,U);
    for (int i = 0; i < dim(); ++i)
    {
        D *= U.at(i,i);
    }
    return D;
}

// REQUIRES: N/A
// MODIFIES: N/A
// EFFECTS:  Returns the transpose of a matrix.
matrix matrix::trans() const
{
    matrix t(column(), row());
    for (int i = 0; i < row(); ++i)
    {
        for (int j = 0; j < column(); ++j)
        {
            t(j,i) = at(i,j);
        }
    }
    return t;
}

// REQUIRES: matrix is square
// MODIFIES: N/A
// EFFECTS:  Returns the adjunct matrix of a given matrix.
matrix matrix::adj() const
{
    matrix A(dim());
    double s;
    for (int i = 0; i < dim(); ++i)
    {
        for (int j = 0; j < dim(); ++j)
        {
            s = ((i + j) % 2 == 0) ? 1: -1;
            A(i,j) = s * cof(i,j).det();
        }
    }
    return A;
}

// REQUIRES: matrix is not singular
// MODIFIES: N/A
// EFFECTS: Returns the inverse matrix of a given matrix;
matrix matrix::inv() const
{
    double D = det();
    if (is_approx(D,0.0)) throw singular_matrix("matrix is singular. Unable "
                                                "to compute inverse.");
    matrix a = adj();
    matrix I(_row);
    for (int i = 0; i < _col; ++i)
    {
        for (int j = 0; j < _row; ++j)
        {
            I(i,j) = adj().at(i,j) / D;
        }
    }
    return I;
}

// REQUIRES: L and U are square matrices with dimension of this matrix
// MODIFIES: The lower triangular matrix is stored in L and the upper triangular
//           matrix is stored in U.
// EFFECTS:  N/A
// NEED TO : Implement generalization for non-square matrices.
//           Implement use of row permutations to catch errors.
void matrix::LU(matrix &L, matrix &U) const
{
    if (!L.is_square() || !U.is_square())
    {
        throw rectangular_matrix("matrix is not square. Unable to factor");
    }
    else if (L.dim() != dim() || U.dim() != dim())
    {
        throw dimension_mismatch("L or U matrix does not have correct "
                                 "dimensions. Unable to factor.");
    }
    
    matrix temp(*this);
    int n = dim();
    matrix C;
    matrix R;
    double pivot;
    
    for (int k = 0; k < n; ++k)
    {
        C = temp.get_col(k);
        R = temp.get_row(k);
        pivot = C(k,0);
        if (!is_approx(pivot,0)) {
            C = C / pivot;
            temp = temp - C*R;
            L.set_col(k,C);
            U.set_row(k,R);
        }
        else
        {
            std::cout << "matrix requires row permutations" << std::endl;
            std::cout << "Step where algorithm failed is " << k << std::endl;
            L = matrix(n);
            U = matrix(n);
            break;
        }
    }
}

// REQUIRES: Q and R are square matrices with dimension of this matrix.
// MODIFIES: The orthonormal matrix is stored in Q and the upper triangular
//           matrix is stored in U.
// EFFECTS:  N/A
void matrix::QR(matrix &Q, matrix &R) const
{
    if (!Q.is_square() || !R.is_square())
    {
        throw rectangular_matrix("matrix is not square. Unable to factor");
    }
    else if (Q.dim() != dim() || R.dim() != dim())
    {
        throw dimension_mismatch("Q or R matrix does not have correct "
                                 "dimensions. Unable to factor.");
    }
    
    matrix Z = zeros(dim(),1);
    matrix sum(Z), u, q, v;
    for (int k = 0; k < column(); ++k)
    {
        u = get_col(k);
        for (int i = 0; i < k; ++i) {
            q = Q.get_col(i);
            sum = (sum + q * (dotProduct(u, q) / dotProduct(q, q)));
        }
        v = u - sum;
        Q.set_col(k, v / norm(v));
        sum = Z;
    }
    R = Q.trans()*(*this);
}

// REQUIRES: row is less than the number of rows in the matrix and greater than
//           0. col is less than the number of columns in the matrix and greater
//           than 0.
// MODIFIES: N/A
// EFFECTS:  Returns the value at the given index.
double matrix::at(const int &row, const int &col) const
{
    if (row < 0 || row > _row-1 || col < 0 || col > _col-1)
    {
        throw invalid_index("Unable to access elements. Indices out of bounds."
                          " row = "+std::to_string(row)+
                          " col = "+std::to_string(col));
    }
    return data[row*column() + col];
}

// REQUIRES: matrix is not empty
// MODIFIES: Deletes dynamically allocated array storing matrix data and
//           sets row, col, and size to 0.
// EFFECTS:  N/A
void matrix::clear()
{
    if (!empty()) delete[] data;
    _row = 0;
    _col = 0;
    _size = 0;
    _empty = true;
    _square = false;
}

// REQUIRES: row and col must be greater than 0
// MODIFIES: Deletes rows or columns if row is less than the number of rows or
//           if col is less than the number of columns. Adds rows or columns
//           initialized to 0 if row is greater than the number of rows or if
//           col is greater than the number of columns.
// EFFECTS:  N/A
void matrix::resize(const int &row, const int &col)
{
    matrix m(row, col);
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            if (i > _row-1 || j > _col-1) m(i,j) = 0;
            else m(i,j) = at(i,j);
        }
    }
    clear();
    _row = m.row();
    _col = m.column();
    _size = m.row()*m.column();
    _square = m.row() == m.column();
    _empty = m.size() > 0;
    data = new double[m.size()];
    for (int i = 0; i < this->row(); ++i)
    {
        for (int j = 0; j < this->column(); ++j)
        {
            (*this)(i,j) = m.at(i,j);
        }
    }
}

// REQUIRES: N/A
// MODIFIES: N/A
// EFFECTS:  Prints matrix to terminal.
void matrix::print(const int &n) const
{
    std::cout.precision(n);
    for (int i = 0; i < row(); ++i) {
        std::cout << "| ";
        for (int j = 0; j < column(); ++j) {
            std::cout << std::setw(n*2) << at(i,j) << ' ';
        }
        std::cout << " |" << std::endl;
    }
    std::cout << std::endl;
}

// REQUIRES: n is greater than or equal to 0 and less than the number of rows
// MODIFIES: N/A
// EFFECTS: Returns the row at the given index
matrix matrix::get_row(const int &n) const
{
    matrix r(1,column());
    for (int i = 0; i < row(); ++i) {
        r(0,i) = at(n,i);
    }
    return r;
}

// REQUIRES: n is greater than or equal to 0 and less than the number of columns
// MODIFIES: N/A
// EFFECTS: Returns the column at the given index
matrix matrix::get_col(const int &n) const
{
    matrix c(row(),1);
    for (int i = 0; i < row(); ++i) {
        c(i,0) = at(i,n);
    }
    return c;
}

// REQUIRES: n is greater than or equal to 0 and less than the number of rows
//           r is a row vector
// MODIFIES: Sets the values of the nth row to the values in r.
// EFFECTS:  N/A
void matrix::set_row(const int &n, const matrix &r)
{
    for (int i = 0; i < row(); ++i) {
        (*this)(n,i) = r.at(0,i);
    }
}

// REQUIRES: n is greater than or equal to 0 and less than the number of columns
//           c is a column vector
// MODIFIES: Sets the values of the nth column to the values in c.
// EFFECTS:  N/A
void matrix::set_col(const int &n, const matrix &c)
{
    for (int i = 0; i < row(); ++i) {
        (*this)(i,n) = c.at(i,0);
    }
}

// REQUIRES: Input matrix must be a row vector
// MODIFIES: Adds row to matrix with given values
// EFFECTS:  N/A
void matrix::push_row_back(const matrix m)
{
    if (empty()) resize(_row+1, m.column());
    else resize(_row+1, _col);
    for (int i = 0; i < _col; ++i) (*this)(_row-1, i) = m.at(0,i);
}

// REQUIRES: Input matrix must be a column vector
// MODIFIES: Adds column to matrix with given values
// EFFECTS:  N/A
void matrix::push_column_back(const matrix m)
{
    if(empty()) resize(m.row(),_col+1);
    else resize(_row, _col+1);
    for (int i = 0; i < _row; ++i) (*this)(i, _col-1) = m.at(i,0);
}

// REQUIRES: There must be at least 1 row in the matrix
// MODIFIES: Deletes bottom row of matrix
// EFFECTS:  N/A
void matrix::pop_row_back()
{
    resize(_row-1, _col);
}

// REQUIRES: There must be at least 1 column in the matrix
// MODIFIES: Deletes bottom column of matrix
// EFFECTS:  N/A
void matrix::pop_column_back()
{
    resize(_row, _col-1);
}

// Non-Member Functions

// REQUIRES: A and B are both column vectors of the same size
// MODIFIES: N/A
// EFFECTS:  Returns the dot product of A and B such that A•B=C
double dotProduct(const matrix &A, const matrix &B)
{
    double C = 0;
    for (int i = 0; i < A.size(); ++i) C += A.at(i,0) * B.at(i,0);
    return C;
}

// REQUIRES: matrix is a column vector
// MODIFIES: N/A
// EFFECTS:  Returns the norm of the vector;
double norm(const matrix &v)
{
    double n = 0;
    for (int i = 0; i < v.size(); ++i) n += v.at(i,0)*v.at(i,0);
    return sqrt(n);
}

// REQUIRES: matrix is a column vector
// MODIFIES: N/A
// EFFECTS:  Returns the norm of the vector;
matrix pow(const matrix &A, const int &n)
{
    matrix m(A);
    double temp;
    for (int i = 0; i < A.row(); ++i)
    {
        for (int j = 0; j < A.column(); ++j)
        {
            if (n == 0) m(i,j) = 1;
            else
            {
                temp = m(i,j);
                for (int k = 0; k < n-1; ++k)
                {
                    m(i,j) *= temp;
                }
            }
        }
    }
    return m;
}

// REQUIRES: U must be an upper triangular matrix. The number of rows in U must
//           be equal to the number of rows in b.
// MODIFIES: N/A
// EFFECTS:  Returns solution to Ux = b;
matrix backwardSub(const matrix &U, const matrix &b)
{
    matrix x(U.column(), b.column());
    int n = U.dim()-1;
    x(n, 0) = b.at(n, 0) / U.at(n, n);
    double sum = 0;
    for (int i = n-1; i >= 0; --i)
    {
        for (int j = n; j > i; --j)
        {
            sum += U.at(i,j)*x.at(j, 0);
        }
        x(i,0) = (b.at(i,0) - sum) / U.at(i,i);
        sum = 0;
    }
    return x;
}

// REQUIRES: L must be an lower triangular matrix. The number of rows in L must
//           be equal to the number of rows in b.
// MODIFIES: N/A
// EFFECTS:  Returns solution to Lx = b;
matrix forwardSub(const matrix &L, const matrix &b)
{
    matrix x(L.column(), b.column());
    double sum = 0;
    x(0, 0) = b.at(0, 0) / L.at(0, 0);
    for (int i = 1; i < L.dim(); ++i)
    {
        for (int j = 0; j < i; ++j)
        {
            sum += L.at(i,j)*x.at(j, 0);
        }
        x(i,0) = (b.at(i,0) - sum) / L.at(i,i);
        sum = 0;
    }
    return x;
}

// REQUIRES: Number of rows in A must be equal to the number of rows in b
// MODIFIES: N/A
// EFFECTS:  Returns solution to the linear equation Ax = b
matrix linearSolve(const matrix &A, const matrix &b,
                      const std::string &type)
{
    double D = A.det();
    if (is_approx(D,0.0))
    {
        throw singular_matrix("matrix is singular. Unable to solve.");
    }
    else if (A.row() != b.row())
    {
        throw std::invalid_argument("Rows in A != rows in b");
    }
    matrix x(A.column(), b.row());
    if (type == "LU")
    {
        matrix L(A.dim()), U(A.dim());
        A.LU(L,U);
        matrix y_star = forwardSub(L,b);
        x = backwardSub(U, y_star);
    }
    else if (type == "QR")
    {
        matrix Q(A.dim()), R(A.dim());
        A.QR(Q,R);
        x = backwardSub(R,Q.trans()*b);
    }
    return x;
}

// REQUIRES: Number of rows in x must be equal to the number of rows in y
// MODIFIES: N/A
// EFFECTS:  Returns coefficient matrix for specified polynomial fit.
// Example:
//        n = 1: Returns coefficients for C1*x + C2 = y
//        n = 2: Returns coefficients for C1*x^2 + C2*x + C3 = y
//        n = m: Returns coefficients for C1*x^m + C2*x^(m-1) + ... + C(m+1) = y
matrix fit(const matrix &x, const matrix &y, const int &n,
              const std::string &type)
{
    matrix phi(x.row(),n+1);
    for (int i = 0; i < n+1; ++i)
    {
        phi.set_col(i,pow(x,n-i));
    }
    matrix A = phi.trans() * phi;
    matrix b = phi.trans() * y;
    return linearSolve(A,b,type);
}

// REQUIRES: order is a combination of X Y and Z. x, y, and z are in radians.
// MODIFIES: N/A
// EFFECTS:  Returns rotation matrix.
matrix rotationmatrix(double x, double y, double z, std::string order)
{
    double _x[9] = {1.0, 0.0, 0.0, 0.0, cos(x), -sin(x), 0.0, sin(x), cos(x)};
    double _y[9] = {cos(y), 0.0, sin(y), 0.0, 1.0, 0.0, -sin(y), 0.0, cos(y)};
    double _z[9] = {cos(z), -sin(z), 0.0, sin(z), cos(z), 0.0, 0.0, 0.0, 1.0};
    
    matrix X(3,3,_x), Y(3,3,_y), Z(3,3,_z);
    if (order == "XYZ") return (X * Y) * Z;
    else if (order == "XZY") return (X * Z) * Y;
    else if (order == "YXZ") return (Y * X) * Z;
    else if (order == "YZX") return (Y * Z) * X;
    else if (order == "ZXY") return (Z * X) * Y;
    else if (order == "ZYX") return (Z * Y) * X;
    else throw std::invalid_argument("Order must be a "
                                     "combination of X Y Z");
}

#endif /* Matrix_h */
