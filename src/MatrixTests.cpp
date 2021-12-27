//
//  matrix_tests.cpp
//  Matrix
//
//  Created by Broderick Riopelle on 12/23/21.
//

#include <stdio.h>
#include "Matrix.h"
#include "unit_test_framework.h"

TEST(A_init_test)
{
    double data[9] = {31, 25, 13, 64, 53, 46, 97, 68, 19};
    Matrix A(3,3,data);
    Matrix B(3,data);
    Matrix C(A);
    Matrix D(3,3);
    Matrix E(3);
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

TEST(B_LU_test)
{
    double data[9] = {31, 25, 13, 64, 53, 46, 97, 68, 19};
    Matrix A(3,data);
    Matrix L(3),U(3);
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

TEST(C_QR_test_1)
{
    double data[9] = {1,1,0,1,2,1,0,3,1};
    Matrix A(3,data);
    Matrix Q(3),R(3);
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

TEST(C_QR_test_2)
{
    double data[16] = {19,9,8,17,16,15,14,3,11,13,12,18,10,7,1,4};
    Matrix A(4,data);
    Matrix Q(4),R(4);
    A.QR(Q,R);
    std::cout << "Q: " << std::endl;
    (Q).print();
    std::cout << "R: " << std::endl;
    R.print();
    std::cout << "Q^T*Q: " << std::endl;
    (Q.trans()*Q).print();
    std::cout << "Q*R: " << std::endl;
    (Q*R).print();
    std::cout << "A: " << std::endl;
    A.print();
}

TEST(D_LinearSolve_test)
{
    double A_data[9] = {31, 25, 13, 64, 53, 46, 97, 68, 19};
    double b_data[3] = {1, 2, 3};
    Matrix A(3,A_data), b(3,1,b_data), x;
    x = linearSolve(A,b);
    x.print();
}

TEST(E_Regression_test_1)
{
    double ans[2] = {2.12280701754,
                     2.33333333333};
    double x_data[5] = {1, 2, 4, 5, 7};
    double y_data[5] = {4, 8, 10, 12, 18};
    Matrix x(5,1,x_data), y(5,1,y_data);
    Matrix C = fit(x,y);
    C.print();
    for(int i = 0; i < 2; ++i) ASSERT_TRUE(is_approx(C.at(i), ans[i]));
}

TEST(F_Regression_test_2)
{
    double ans[3] = {1.35780885781,
                    -4.60536130536,
                     19.4125874126};
    double x_data[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    double y_data[11] = {20, 15, 16, 18, 23, 30, 40, 55, 70, 86, 110};
    Matrix x(11,1,x_data), y(11,1,y_data);
    Matrix C = fit(x,y,2);
    C.print();
    for(int i = 0; i < 3; ++i) ASSERT_TRUE(is_approx(C.at(i), ans[i]));
}

TEST(G_Regression_test_3)
{
    double ans[4] = {-0.262820512821,
                      4.71503496503,
                     -16.9696969697,
                      32.9160839161};
    double x_data[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    double y_data[11] = {33, 20, 16, 18, 23, 33, 44, 55, 65, 70, 72};
    Matrix x(11,1,x_data), y(11,1,y_data);
    Matrix C = fit(x,y,3);
    C.print();
    for(int i = 0; i < 3; ++i) ASSERT_TRUE(is_approx(C.at(i), ans[i]));
}

TEST(H_Rotation_test)
{
    double ans[9] = {0.2919266, -0.4546487,  0.8414710,
                     0.8372224, -0.3038967, -0.4546487,
                     0.4624257,  0.8372224,  0.2919266};
    Matrix R = rotationMatrix(1,1,1);
    R.print();
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            ASSERT_TRUE(is_approx(R.at(i,j), ans[j + 3*i], 1e-6));
        }
    }
}

TEST_MAIN()
