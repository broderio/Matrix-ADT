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

class matrix {
private:
    int _row;       // number of rows
    int _col;       // number of columns
    int _size;      // number of elements in matrix
    bool _square;   // boolean for square matrix
    bool _empty;    // boolean for empty matrix
    double* data;   // dynamically allocated array with data
public:
    // Constructors
    matrix();
    matrix(const int &dim);
    matrix(const int &dim, const double *arr);
    matrix(const int &row_size, const int &col_size);
    matrix(const int &row_size, const int &col_size, const double *arr);
    
    // Big Three
    matrix(const matrix &other);
    ~matrix();
    const matrix& operator=(const matrix &rhs);
    
    // Overloaded Operators
    double &operator()(const int &row, const int &col); // Element Access
    matrix operator+(const matrix &rhs) const;          // matrix Addition
    matrix operator+(const double &a) const;            // Scalar Addition
    matrix operator-(const matrix &rhs) const;          // matrix Subtraction
    matrix operator-(const double &a) const;            // Scalar Subtraction
    matrix operator*(const matrix &rhs) const;          // matrix Multiplication
    matrix operator*(const double &a) const;            // Scalar Multiplication
    matrix operator/(const double &a) const;            // Scalar Division
    
    // Functions
    int size() const;
    int row() const;
    int column() const;
    
    bool is_square() const;
    bool empty() const;
    int dim() const;
    bool contains(const double &n) const;
    
    matrix cof(const int &row, const int &col) const; // Cofactor matrix
    double det() const;                               // Determinant
    matrix trans() const;                             // Transpose
    matrix adj() const;                               // Adjugate
    matrix inv() const;                               // Inverse
    
    void LU(matrix &L, matrix &U) const;              // LU Decomposition
    void QR(matrix &L, matrix &U) const;              // QR Decomposition
    
    double at(const int &row, const int &col = 0) const;
    void clear();
    void resize(const int &row, const int &col);
    
    matrix get_row(const int &n) const;
    matrix get_col(const int &n) const;
    void set_row(const int &n, const matrix &r);
    void set_col(const int &n, const matrix &c);
    
    void push_row_back(const matrix m);
    void push_column_back(const matrix m);
    void pop_row_back();
    void pop_column_back();
    
    void print(const int &n = 5) const;
};

double dotProduct(const matrix &A, const matrix &B);
double norm(const matrix &v);
matrix pow(const matrix &A, const int &n);

matrix backwardSub(const matrix &U, const matrix &b);
matrix forwardSub(const matrix &L, const matrix &b);

matrix linearSolve(const matrix &A, const matrix &b,
                   const std::string &type = "LU");
matrix fit(const matrix &x, const matrix &y, const int &n=1,
           const std::string &type = "LU");

matrix rotationmatrix(double x, double y, double z, std::string order="XYZ");

matrix zeros(const int &dim);
matrix zeros(const int &row, const int &col);
matrix ones(const int &dim);
matrix ones(const int &row, const int &col);

#endif /*MatrixClassDef_h */
