//
//  MatrixClassDef.h
//  Matrix
//
//  Created by Broderick Riopelle on 12/26/21.
//

#ifndef MatrixClassDef_h
#define MatrixClassDef_h

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

#endif /* MatrixClassDef_h */
