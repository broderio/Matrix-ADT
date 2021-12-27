//
//  MatrixClassDef.h
//  Matrix
//
//  Created by Broderick Riopelle on 12/26/21.
//

#ifndef MatrixClassDef_h
#define MatrixClassDef_h

const double EPSILON = 1e-10;
bool is_approx(const double &a, const double &b, const double &tol = EPSILON);

class Matrix {
private:
    int _row; // number of rows
    int _col; // number of columns
    int _size; // number of elements in matrix
    bool _square;
    bool _empty;
    double* data; // dynamically allocated array with data
public:
    // Constructors
    Matrix();
    Matrix(const int &dim);
    Matrix(const int &dim, const double *arr);
    Matrix(const int &row_size, const int &col_size);
    Matrix(const int &row_size, const int &col_size, const double *arr);
    
    // Big Three
    Matrix(const Matrix &other);
    ~Matrix();
    const Matrix& operator=(const Matrix &rhs);
    
    // Overloaded Operator
    double &operator()(const int &row, const int &col);
    Matrix operator+(const Matrix &rhs) const;
    Matrix operator+(const double &a) const;
    Matrix operator-(const Matrix &rhs) const;
    Matrix operator*(const Matrix &rhs) const;
    Matrix operator*(const double &a) const;
    Matrix operator/(const double &a) const;
    
    // Functions
    int size() const;
    int row() const;
    int column() const;
    
    bool is_square() const;
    bool empty() const;
    int dim() const;
    
    Matrix cof(const int &row, const int &col) const; // Cofactor Matrix
    double det() const;                                       // Determinant
    Matrix trans() const;                             // Transpose
    Matrix adj() const;                               // Adjugate
    Matrix inv() const;                           // Inverse
    
    void LU(Matrix &L, Matrix &U) const;           // LU Decomposition
    void QR(Matrix &L, Matrix &U) const;           // QR Decomposition
    
    double at(const int &row, const int &col = 0) const;
    void clear();
    void resize(const int &row, const int &col);
    
    Matrix get_row(const int &n) const;
    Matrix get_col(const int &n) const;
    void set_row(const int &n, const Matrix &r);
    void set_col(const int &n, const Matrix &c);
    
    void push_row_back(const Matrix m);
    void push_column_back(const Matrix m);
    void pop_row_back();
    void pop_column_back();
    
    void print(const int &n = 5) const;
};

double dotProduct(const Matrix &A, const Matrix &B);
double norm(const Matrix &v);
Matrix pow(const Matrix &A, const int &n);

Matrix backwardSub(const Matrix &U, const Matrix &b);
Matrix forwardSub(const Matrix &L, const Matrix &b);

Matrix linearSolve(const Matrix &A, const Matrix &b, const std::string &type = "LU");
Matrix fit(const Matrix &x, const Matrix &y, const int &n=1, const std::string &type = "LU");

Matrix rotationMatrix(double x, double y, double z, std::string order);

Matrix zeros(const int &dim);
Matrix zeros(const int &row, const int &col);
Matrix ones(const int &dim);
Matrix ones(const int &row, const int &col);

#endif /* MatrixClassDef_h */
