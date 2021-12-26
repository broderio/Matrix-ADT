//
//  matrix_tests.cpp
//  Matrix
//
//  Created by Broderick Riopelle on 12/23/21.
//

#include <stdio.h>
#include "Matrix.hpp"
#include "unit_test_framework.h"

TEST(init_test)
{
    float data[9] = {31, 25, 13, 64, 53, 46, 97, 68, 19};
    Matrix<float> A(3,3,data);
    Matrix<float> B(3,data);
    Matrix<float> C(A);
    Matrix<float> D(3,3);
    Matrix<float> E(3);
    std::cout << "A: " << std::endl;
    A.print();
    std::cout << "B: " << std::endl;
    B.print();
    std::cout << "C: " << std::endl;
    C.print();
    ASSERT_FALSE(A.empty());
    ASSERT_FALSE(B.empty());
    ASSERT_FALSE(C.empty());
    ASSERT_FALSE(D.empty());
    ASSERT_FALSE(E.empty());
    ASSERT_TRUE(A.is_square());
    ASSERT_TRUE(B.is_square());
    ASSERT_TRUE(C.is_square());
    ASSERT_TRUE(D.is_square());
    ASSERT_TRUE(E.is_square());
    ASSERT_EQUAL(A.size(), 9);
    ASSERT_EQUAL(B.size(), 9);
    ASSERT_EQUAL(C.size(), 9);
    ASSERT_EQUAL(D.size(), 9);
    ASSERT_EQUAL(E.size(), 9);
    ASSERT_EQUAL(A.row(), 3);
    ASSERT_EQUAL(B.row(), 3);
    ASSERT_EQUAL(C.row(), 3);
    ASSERT_EQUAL(D.row(), 3);
    ASSERT_EQUAL(E.row(), 3);
    ASSERT_EQUAL(A.column(), 3);
    ASSERT_EQUAL(B.column(), 3);
    ASSERT_EQUAL(C.column(), 3);
    ASSERT_EQUAL(D.column(), 3);
    ASSERT_EQUAL(E.column(), 3);
}

TEST(LU_test)
{
    float data[9] = {31, 25, 13, 64, 53, 46, 97, 68, 19};
    Matrix<float> A(3,data);
    Matrix<float> L(3),U(3);
    A.LU(L,U);
    std::cout << "L: " << std::endl;
    L.print();
    std::cout << "U: " << std::endl;
    U.print();
    std::cout << "L*U: " << std::endl;
    (L*U).print();
    std::cout << "A: " << std::endl;
    A.print();
}

TEST(QR_test)
{
    float data[9] = {1,1,0,1,2,1,0,3,1};
    Matrix<float> A(3,data);
    Matrix<float> Q(3),R(3);
    A.QR(Q,R);
    std::cout << "Q: " << std::endl;
    (Q).print();
    std::cout << "R: " << std::endl;
    R.print();
    std::cout << "Q*R: " << std::endl;
    (Q*R).print();
    std::cout << "A: " << std::endl;
    A.print();
}

TEST(Regression_test)
{
    float x_data[5] = {1, 2, 4, 5, 7};
    float y_data[5] = {4, 8, 10, 12, 18};
    Matrix<float> x(5,1,x_data), y(5,1,y_data);
}

TEST(LinearSolve_test)
{
    float A_data[9] = {31, 25, 13, 64, 53, 46, 97, 68, 19};
    float b_data[3] = {1, 2, 3};
    Matrix<float> A(3,A_data), b(3,1,b_data), x;
    x = linearSolve(A,b);
    x.print();
}

TEST_MAIN()
