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

template <typename T>
class Matrix {
private:
    int _row; // number of rows
    int _col; // number of columns
    int _size; // number of elements in matrix
    bool _square;
    bool _empty;
    T* data; // dynamically allocated array with data
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
    Matrix<double> inv() const;                           // Inverse
    
    void LU(Matrix<T> &L, Matrix<T> &U) const;           // LU Decomposition
    void QR(Matrix<T> &L, Matrix<T> &U) const;           // QR Decomposition
    
    T at(const int &row, const int &col = 0) const;
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

template <typename T>
T dotProduct(const Matrix<T> &A, const Matrix<T> &B);

template <typename T>
double norm(const Matrix<T> &v);

template <typename T>
Matrix<T> pow(const Matrix<T> &A, const int &n);

template <typename T>
Matrix<T> backwardSub(const Matrix<T> &U, const Matrix<T> &b);

template <typename T>
Matrix<T> forwardSub(const Matrix<T> &L, const Matrix<T> &b);

template <typename T>
Matrix<T> linearSolve(const Matrix<T> &A, const Matrix<T> &b, const std::string &type = "LU");

template <typename T>
Matrix<T> fit(const Matrix<T> &x, const Matrix<T> &y, const int &n=1, const std::string &type = "LU");

Matrix<double> rotationMatrix(double x, double y, double z, std::string order);

Matrix<double> zeros(const int &dim);

Matrix<double> zeros(const int &row, const int &col);

Matrix<double> ones(const int &dim);

Matrix<double> ones(const int &row, const int &col);

#endif /* MatrixClassDef_h */
