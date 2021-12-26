//
//  Matrix.hpp
//  Matrix
//
//  Created by Broderick Riopelle on 12/22/21.
//

#ifndef Matrix_hpp
#define Matrix_hpp

#include <stdio.h>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class dimension_error {
private:
    std::string msg;
public:
    dimension_error(std::string msg_in);
    std::string what() const;
};

class index_error {
private:
    std::string msg;
public:
    index_error(std::string msg_in);
    std::string what() const;
};

class not_square_error {
private:
    std::string msg;
public:
    not_square_error(std::string msg_in);
    std::string what() const;
};

class is_singular_error {
private:
    std::string msg;
public:
    is_singular_error(std::string msg_in);
    std::string what() const;
};

class is_empty_error {
private:
    std::string msg;
public:
    is_empty_error(std::string msg_in);
    std::string what() const;
};

template <typename T>
class Matrix {
private:
    int _row; // number of rows
    int _col; // number of columns
    int _size; // number of elements in matrix
    bool _square;
    bool _empty;
    T* data; // dynamically allocated array with data
    const float EPSILON = 1e-8;
    bool float_equal(const float &a, const float &b) const;
public:
    // Constructors
    Matrix();
    Matrix(const int &dim);
    Matrix(const int &dim, const T *arr);
    Matrix(const int &row_size, const int &col_size);
    Matrix(const int &row_size, const int &col_size, const T *arr);
    
    // Big Three
    Matrix(const Matrix<T> &other);
    ~Matrix();
    const Matrix<T>& operator=(const Matrix<T> &rhs);
    
    // Overloaded Operator
    T &operator()(const int &row, const int &col);
    Matrix<T> operator+(const Matrix<T> &rhs) const;
    Matrix<T> operator+(const T &a) const;
    Matrix<T> operator-(const Matrix<T> &rhs) const;
    Matrix<T> operator*(const Matrix<T> &rhs) const;
    Matrix<T> operator*(const T &a) const;
    Matrix<T> operator/(const T &a) const;
    
    // Functions
    int size() const;
    int row() const;
    int column() const;
    
    bool is_square() const;
    bool empty() const;
    int dim() const;
    
    Matrix<T> cof(const int &row, const int &col) const; // Cofactor Matrix
    T det() const;                                       // Determinant
    Matrix<T> trans() const;                             // Transpose
    Matrix<T> adj() const;                               // Adjugate
    Matrix<float> inv() const;                           // Inverse
    
    void LU(Matrix<T> &L, Matrix<T> &U) const;           // LU Decomposition
    void QR(Matrix<T> &L, Matrix<T> &U) const;           // QR Decomposition
    
    T at(const int &row, const int &col) const;
    void clear();
    void resize(const int &row, const int &col);
    
    Matrix<T> get_row(const int &n) const;
    Matrix<T> get_col(const int &n) const;
    void set_row(const int &n, const Matrix<T> &r);
    void set_col(const int &n, const Matrix<T> &c);
    
    void push_row_back(const Matrix<T> m);
    void push_column_back(const Matrix<T> m);
    void pop_row_back();
    void pop_column_back();
    
    void print(const int &n = 5) const;
};

// Error Messages

dimension_error::dimension_error(std::string msg_in)
: msg(msg_in) {}

std::string dimension_error::what() const
{
    return msg;
}

index_error::index_error(std::string msg_in)
: msg(msg_in) {}

std::string index_error::what() const
{
    return msg;
}

not_square_error::not_square_error(std::string msg_in)
: msg(msg_in) {}

std::string not_square_error::what() const
{
    return msg;
}

is_singular_error::is_singular_error(std::string msg_in)
: msg(msg_in) {}

std::string is_singular_error::what() const
{
    return msg;
}

is_empty_error::is_empty_error(std::string msg_in)
: msg(msg_in) {}

std::string is_empty_error::what() const
{
    return msg;
}

// REQUIRES: N/A
// MODIFIES: N/A
// EFFECTS: Returns true if a is within tolerance of b.
template <typename T>
bool Matrix<T>::float_equal(const float &a, const float &b) const
{
    if (a-b <= EPSILON) return true;
    return false;
}

// Constructors

// REQUIRES: N/A
// MODIFIES: Initializes 0 by 0 empty matrix
// EFFECTS: N/A
template <typename T>
Matrix<T>::Matrix()
: _row(0), _col(0), _size(0), data(nullptr), _square(false), _empty(true) {}

// REQUIRES: N/A
// MODIFIES: Initializes empty square matrix with dimension dim
// EFFECTS: N/A
template <typename T>
Matrix<T>::Matrix(const int &dim)
: _row(dim), _col(dim), _size(dim*dim), data(new T[dim*dim]), _square(true),
_empty(false) {}

// REQUIRES: arr must have dim^2 elements
// MODIFIES: Initializes square matrix with dimension dim and with data in arr.
//           The data is assumed to be in order of a 2D matrix starting left to
//           right and going downward.
//           Example: Matrix(2, {1, 2, 3, 4}) would result in [ 1 2 ]
//                                                            [ 3 4 ]
// EFFECTS: N/A
template <typename T>
Matrix<T>::Matrix(const int &dim, const T *arr)
: _row(dim), _col(dim), _size(dim*dim), data(new T[dim*dim]), _square(true),
_empty(false)
{
    for (int i = 0; i < size(); ++i) data[i] = arr[i];
}

// REQUIRES: N/A
// MODIFIES: Initializes row_size by col_size empty matrix
// EFFECTS: N/A
template <typename T>
Matrix<T>::Matrix(const int &row_size, const int &col_size)
: _row(row_size), _col(col_size), _size(row_size*col_size),
data(new T[row_size*col_size]), _square(_row == _col), _empty(false)
{}

// REQUIRES: arr must be of size row_size * col_size
// MODIFIES: Initializes row_size by col_size matrix with data in arr. The data
//           is assumed to be in order of a 2D matrix starting left to right and
//           going downward.
//           Example: Matrix(2, 4, {1, 2, 3, 4, 5, 6, 7, 8}) would result in
//                    [ 1 2 3 4 ]
//                    [ 5 6 7 8 ]
// EFFECTS: N/A
template <typename T>
Matrix<T>::Matrix(const int &row_size, const int &col_size, const T *arr)
: _row(row_size), _col(col_size), _size(row_size*col_size),
  data(new T[row_size*col_size]), _square(_row == _col), _empty(false)
{
      for (int i = 0; i < size(); ++i) data[i] = arr[i];
}


// Big Three

// REQUIRES: N/A
// MODIFIES: Creates a deep copy of other.
// EFFECTS:  N/A
template <typename T>
Matrix<T>::Matrix(const Matrix<T> &other)
: _row(other.row()), _col(other.column()), _size(other.size()), _empty(other.empty()),
  data(new T[other.size()]), _square(other.is_square())
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
template <typename T>
Matrix<T>::~Matrix()
{
    clear();
}

// REQUIRES: N/A
// MODIFIES: Deletes original matrix and creates a deep copy of rhs.
// EFFECTS:  Returns the matrix after operation.
template <typename T>
const Matrix<T> &Matrix<T>::operator=(const Matrix<T> &rhs)
{
    if (this->data == rhs.data) return *this;
    clear();
    _row = rhs.row();
    _col = rhs.column();
    _size = rhs.size();
    data = new T[_size];
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
template <typename T>
T &Matrix<T>::operator()(const int &row, const int &col)
{
    if (row < 0 || row > _row-1 || col < 0 || col > _col-1)
    {
        throw index_error("Unable to access elements. Indices out of bounds."
                          " row = "+std::to_string(row)+
                          " col = "+std::to_string(col));
    }
    return data[row*column() + col];
}

// REQUIRES: rhs is the same dimensions as lhs
// MODIFIES: N/A
// EFFECTS:  Returns the sum of two matrices
template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &rhs) const
{
    if (rhs.column() != column() || rhs.row() != row())
    {
        throw dimension_error("Dimension mismatch. Matrix addition"
                              "requires matrices of equal dimensions");
    }
    Matrix<T> m(*this);
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
template <typename T>
Matrix<T> Matrix<T>::operator+(const T &a) const
{
    Matrix<T> m(*this);
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
template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &rhs) const
{
    if (rhs.column() != column() || rhs.row() != row())
    {
        throw dimension_error("Dimensions mismatch. Matrix subtraction"
                              "requires matrices of equal dimensions");
    }
    Matrix<T> m(*this);
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
template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &rhs) const
{
    if (rhs.row() != column())
    {
        throw dimension_error("Dimensions mismatch. Matrix multiplication "
                              "requires number of columns in first matrix "
                              "to equal number of rows in second matrix.");
    }
    Matrix<T> m(_row, rhs.column());
    for (int i = 0; i < _row; ++i)
    {
        for (int j = 0; j < rhs.column(); ++j)
        {
            T sum = 0;
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
template <typename T>
Matrix<T> Matrix<T>::operator*(const T &a) const
{
    Matrix<T> m(*this);
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
template <typename T>
Matrix<T> Matrix<T>::operator/(const T &a) const
{
    Matrix<T> m(*this);
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
template <typename T>
int Matrix<T>::size() const
{
    return _size;
}

// REQUIRES: N/A
// MODIFIES: N/A
// EFFECTS:  Returns number of rows in matrix
template <typename T>
int Matrix<T>::row() const
{
    return _row;
}

// REQUIRES: N/A
// MODIFIES: N/A
// EFFECTS:  Returns number of columns in matrix
template <typename T>
int Matrix<T>::column() const
{
    return _col;
}

// REQUIRES: N/A
// MODIFIES: N/A
// EFFECTS:  Returns true if matrix is square, false otherwise
template <typename T>
bool Matrix<T>::is_square() const
{
    return _square;
}

// REQUIRES: N/A
// MODIFIES: N/A
// EFFECTS:  Returns true if matrix is empty, false otherwise
template <typename T>
bool Matrix<T>::empty() const
{
    return _empty;
}

// REQUIRES: N/A
// MODIFIES: N/A
// EFFECTS:  Returns dimension if matrix is square, 0 otherwise
template <typename T>
int Matrix<T>::dim() const
{
    if (_square) return _row;
    else return 0;
}

// REQUIRES: Matrix must be square.
// MODIFIES: N/A
// EFFECTS:  Returns cofactor matrix of given index.
template <typename T>
Matrix<T> Matrix<T>::cof(const int &row, const int &col) const
{
    if (!is_square()) throw not_square_error("Matrix is not square. Unable to "
                                             "to compute cofactor.");
    Matrix<T> c(this->row()-1);
    for (int i = 0; i < c.row(); ++i)
    {
        for (int j = 0; j < c.column(); ++j)
        {
            if (i != row && j != col) c(i,j) = at(i,j);
        }
    }
    return c;
}

// REQUIRES: Matrix must be square
// MODIFIES: N/A
// EFFECTS:  Returns determinant of matrix
template <typename T>
T Matrix<T>::det() const
{
    if (!is_square()) throw not_square_error("Matrix is not square. Unable to "
                                             "to compute determinant.");
    Matrix<T> L,U;
    float D = 1;
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
template <typename T>
Matrix<T> Matrix<T>::trans() const
{
    Matrix<T> t(column(), row());
    for (int i = 0; i < row(); ++i)
    {
        for (int j = 0; j < column(); ++j)
        {
            t(j,i) = at(i,j);
        }
    }
    return t;
}

// REQUIRES: Matrix is square
// MODIFIES: N/A
// EFFECTS:  Returns the adjunct matrix of a given matrix.
template <typename T>
Matrix<T> Matrix<T>::adj() const
{
    Matrix<T> A(dim());
    T s;
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

// REQUIRES: Matrix is not singular
// MODIFIES: N/A
// EFFECTS: Returns the inverse matrix of a given matrix;
template <typename T>
Matrix<float> Matrix<T>::inv() const
{
    float D = det();
    if (float_equal(D,0.0)) throw is_singular_error("Matrix is singular. Unable "
                                                "to compute inverse.");
    Matrix<T> a = adj();
    Matrix<float> I(_row);
    for (int i = 0; i < _col; ++i)
    {
        for (int j = 0; j < _row; ++j)
        {
            I(i,j) = adj(i,j) / D;
        }
    }
}


// REQUIRES: L and U are square matrices with dimension of this matrix
// MODIFIES: The lower triangular matrix is stored in L and the upper triangular
//           matrix is stored in U.
// EFFECTS:  N/A
template <typename T>
void Matrix<T>::LU(Matrix<T> &L, Matrix<T> &U) const
{
    Matrix<T> temp(*this);
    int n = dim();
    Matrix<T> C;
    Matrix<T> R;
    float pivot;
    
    for (int k = 0; k < n; ++k)
    {
        C = temp.get_col(k);
        R = temp.get_row(k);
        pivot = C(k,0);
        if (!float_equal(pivot,0)) {
            C = C / pivot;
            temp = temp - C*R;
            L.set_col(k,C);
            U.set_row(k,R);
        }
        else
        {
            std::cout << "Matrix requires row permutations" << std::endl;
            std::cout << "Step where algorithm failed is " << k << std::endl;
            L = Matrix<T>(n);
            U = Matrix<T>(n);
            break;
        }
    }
}

// REQUIRES: Q and R are square matrices with dimension of this matrix.
// MODIFIES: The orthonormal matrix is stored in Q and the upper triangular
//           matrix is stored in U.
// EFFECTS:  N/A
template <typename T>
void Matrix<T>::QR(Matrix<T> &Q, Matrix<T> &R) const
{
    float init[3] = {0, 0, 0};
    Matrix<T> zeros(3,1,init), sum(zeros), u, q, v;
    for (int k = 0; k < column(); ++k)
    {
        u = get_col(k);
        for (int i = 0; i < k; ++i) {
            q = Q.get_col(i);
            sum = (sum + q * (dotProduct(u, q) / dotProduct(q, q)));
        }
        v = u - sum;
        Q.set_col(k, v / norm(v));
        sum = zeros;
    }
    R = Q.trans()*(*this);
}

// REQUIRES: row is less than the number of rows in the matrix and greater than
//           0. col is less than the number of columns in the matrix and greater
//           than 0.
// MODIFIES: N/A
// EFFECTS:  Returns the value at the given index.
template <typename T>
T Matrix<T>::at(const int &row, const int &col) const
{
    if (row < 0 || row > _row-1 || col < 0 || col > _col-1)
    {
        throw index_error("Unable to access elements. Indices out of bounds."
                          " row = "+std::to_string(row)+
                          " col = "+std::to_string(col));
    }
    return data[row*column() + col];
}

// REQUIRES: Matrix is not empty
// MODIFIES: Deletes dynamically allocated array storing matrix data and
//           sets row, col, and size to 0.
// EFFECTS:  N/A
template <typename T>
void Matrix<T>::clear()
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
template <typename T>
void Matrix<T>::resize(const int &row, const int &col)
{
    Matrix<T> m(row, col);
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
    data = new T[m.size()];
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
template <typename T>
void Matrix<T>::print(const int &n) const
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
template <typename T>
Matrix<T> Matrix<T>::get_row(const int &n) const
{
    Matrix<T> r(1,column());
    for (int i = 0; i < row(); ++i) {
        r(0,i) = at(n,i);
    }
    return r;
}

// REQUIRES: n is greater than or equal to 0 and less than the number of columns
// MODIFIES: N/A
// EFFECTS: Returns the column at the given index
template <typename T>
Matrix<T> Matrix<T>::get_col(const int &n) const
{
    Matrix<T> c(row(),1);
    for (int i = 0; i < row(); ++i) {
        c(i,0) = at(i,n);
    }
    return c;
}

// REQUIRES: n is greater than or equal to 0 and less than the number of rows
//           r is a row vector
// MODIFIES: Sets the values of the nth row to the values in r.
// EFFECTS:  N/A
template <typename T>
void Matrix<T>::set_row(const int &n, const Matrix<T> &r)
{
    for (int i = 0; i < row(); ++i) {
        (*this)(n,i) = r.at(0,i);
    }
}

// REQUIRES: n is greater than or equal to 0 and less than the number of columns
//           c is a column vector
// MODIFIES: Sets the values of the nth column to the values in c.
// EFFECTS:  N/A
template <typename T>
void Matrix<T>::set_col(const int &n, const Matrix<T> &c)
{
    for (int i = 0; i < row(); ++i) {
        (*this)(i,n) = c.at(i,0);
    }
}

// REQUIRES: Input matrix must be a row vector
// MODIFIES: Adds row to matrix with given values
// EFFECTS:  N/A
template <typename T>
void Matrix<T>::push_row_back(const Matrix<T> m)
{
    if (empty()) resize(_row+1, m.column());
    else resize(_row+1, _col);
    for (int i = 0; i < _col; ++i) (*this)(_row-1, i) = m.at(0,i);
}

// REQUIRES: Input matrix must be a column vector
// MODIFIES: Adds column to matrix with given values
// EFFECTS:  N/A
template <typename T>
void Matrix<T>::push_column_back(const Matrix<T> m)
{
    if(empty()) resize(m.row(),_col+1);
    else resize(_row, _col+1);
    for (int i = 0; i < _row; ++i) (*this)(i, _col-1) = m.at(i,0);
}

// REQUIRES: There must be at least 1 row in the matrix
// MODIFIES: Deletes bottom row of matrix
// EFFECTS:  N/A
template <typename T>
void Matrix<T>::pop_row_back()
{
    resize(_row-1, _col);
}

// REQUIRES: There must be at least 1 column in the matrix
// MODIFIES: Deletes bottom column of matrix
// EFFECTS:  N/A
template <typename T>
void Matrix<T>::pop_column_back()
{
    resize(_row, _col-1);
}

// Non-Member Functions

// REQUIRES: A and B are both column vectors of the same size
// MODIFIES: N/A
// EFFECTS:  Returns the dot product of A and B such that A•B=C
template <typename T>
T dotProduct(const Matrix<T> &A, const Matrix<T> &B)
{
    T C = 0;
    for (int i = 0; i < A.size(); ++i) C += A.at(i,0) * B.at(i,0);
    return C;
}

// REQUIRES: Matrix is a column vector
// MODIFIES: N/A
// EFFECTS:  Returns the norm of the vector;
template <typename T>
float norm(const Matrix<T> &v)
{
    float n = 0;
    for (int i = 0; i < v.size(); ++i) n += v.at(i,0)*v.at(i,0);
    return sqrt(n);
}

// REQUIRES: U must be an upper triangular matrix. The number of rows in U must
//           be equal to the number of rows in b.
// MODIFIES: N/A
// EFFECTS:  Returns solution to Ux = b;
template <typename T>
Matrix<T> backwardSub(const Matrix<T> &U, const Matrix<T> &b)
{
    Matrix<T> x(U.column(), b.column());
    int n = U.dim()-1;
    x(n, 0) = b.at(n, 0) / U.at(n, n);
    for (int i = n-1; i >= 0; --i)
    {
        x(i,0) = (b.at(i,0) - U.at(i,i+1)*x.at(i+1, 0)) / U.at(i,i);
    }
    return x;
}

// REQUIRES: L must be an lower triangular matrix. The number of rows in L must
//           be equal to the number of rows in b.
// MODIFIES: N/A
// EFFECTS:  Returns solution to Lx = b;
template <typename T>
Matrix<T> forwardSub(const Matrix<T> &L, const Matrix<T> &b)
{
    Matrix<T> x(L.column(), b.column());
    x(0, 0) = b.at(0, 0) / L.at(0, 0);
    for (int i = 1; i < L.dim(); ++i)
    {
        x(i,0) = (b.at(i,0) - L.at(i,i+1)*x.at(i-1, 0)) / L.at(i,i);
    }
    return x;
}

// REQUIRES: Number of rows in A must be equal to the number of rows in b
// MODIFIES: N/A
// EFFECTS:  Returns solution to the linear equation Ax = b
template <typename T>
Matrix<T> linearSolve(const Matrix<T> &A, const Matrix<T> &b)
{
    Matrix<T> x(A.column(), b.row()), Q(A.dim()), R(A.dim());
    (A.trans()*A).QR(Q,R);
    x = backwardSub(R,Q.trans()*(A.trans()*b));
    return x;
}

// REQUIRES: order is a combination of X Y and Z. x, y, and z are in radians.
// MODIFIES: N/A
// EFFECTS:  Returns rotation matrix.
Matrix<float> rotationMatrix(float x, float y, float z, std::string order)
{
    float _x[9] = {1.0, 0.0, 0.0, 0.0, cos(x), -sin(x), 0.0, sin(x), cos(x)};
    float _y[9] = {cos(y), 0.0, sin(y), 0.0, 1.0, 0.0, -sin(y), 0.0, cos(y)};
    float _z[9] = {cos(z), -sin(z), 0.0, sin(z), cos(z), 0.0, 0.0, 0.0, 1.0};
    
    Matrix<float> X(3,3,_x), Y(3,3,_y), Z(3,3,_z);
    if (order == "XYZ") return X * Y * Z;
    else if (order == "XZY") return X * Z * Y;
    else if (order == "YXZ") return Y * X * Z;
    else if (order == "YZX") return Y * Z * X;
    else if (order == "ZXY") return Z * X * Y;
    else if (order == "ZYX") return Z * Y * X;
    else throw std::invalid_argument("Order must be a "
                                     "combination of X Y Z");
}

#endif /* Matrix_hpp */