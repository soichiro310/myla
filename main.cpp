//
//  main.cpp
//  myla
//
//  Created by 佐藤荘一朗 on 2018/07/12.
//  Copyright © 2018年 佐藤荘一朗. All rights reserved.
//

#include "myla.h"
#include <ctime>
#include <random>

#define N 10

void kadai1(){
    clock_t tik, tok;
    Matrix3d A;
    A <<
    2, 3, 0,
    0, 2, 1,
    5, -4, -4;
    
    cout << "行列A" << endl;
    cout << A << endl;
    tik = clock();
    cout << "行列式(Mydet)：" << MyDet(A) << endl;
    tok = clock();
    cout << "実行時間：" << (double)(tok - tik)/CLOCKS_PER_SEC << "[s]"  << endl;
    
    cout << "--------------------------------" <<  endl;
    
    tik = clock();
    cout << "行列式(Eigen)：" << A.determinant() << endl;
    tok = clock();
    cout << "実行時間：" << (double)(tok - tik)/CLOCKS_PER_SEC << "[s]"  << endl;
    
}

void kadai2(){
    clock_t tik, tok;
    Matrix3d A;
    Vector3d b;
    A <<
    3, -4, 5,
    -7, 8, -9,
    11, -5, 6;
    
    b <<
    4,
    -4,
    -3;
    
    cout << "行列A" << endl;
    cout << A << endl;
    cout << "解ベクトルb" << endl;
    cout << b << endl;
    tik = clock();
    cout << "方程式の解(MyCramer) x = " << endl << MyCramer(A, b) << endl;
    tok = clock();
    cout << "実行時間：" << (double)(tok - tik)/CLOCKS_PER_SEC << "[s]"  << endl;
    
    cout << "--------------------------------" <<  endl;
    
    tik = clock();
    cout << "方程式の解(Eigen) x = " << endl << A.fullPivLu().solve(b) << endl;
    tok = clock();
    cout << "実行時間：" << (double)(tok - tik)/CLOCKS_PER_SEC << "[s]"  << endl;
    
}

void kadai3(){
    clock_t tik, tok;
    Matrix3d A;
    Vector3d eigValue;
    Matrix3d eigVector;
    
    A <<
    1, 0, -1,
    1, 2, 1,
    2, 2, 3;
    
    cout << "行列A" << endl;
    cout << A << endl;
    tik = clock();
    eigValue = CalcEigenValue(A);
    eigVector = CalcEigenVector(A, eigValue);
    cout << "固有値(MyFunc)" << endl << eigValue << endl;
    cout << "固有ベクトル(Myfunc)" << endl << eigVector << endl;
    tok = clock();
    cout << "実行時間：" << (double)(tok - tik)/CLOCKS_PER_SEC << "[s]"  << endl;
    
    cout << "--------------------------------" <<  endl;
    
    tik = clock();
    EigenSolver<Matrix<double,3,3> > es(A);
    
    cout << "固有値(Eigen)" << endl << es.eigenvalues() << endl;
    cout << "固有ベクトル(Eigen)" << endl << es.eigenvectors() << endl;
    tok = clock();
    cout << "実行時間：" << (double)(tok - tik)/CLOCKS_PER_SEC << "[s]" << endl;
}

void kadai4(){
    clock_t tik, tok;
    MatrixXd A, eigA;
    VectorXd b;
    VectorXd eigValue;
    MatrixXd eigVector;
    A.resize(N, N);
    eigA.resize(N, N);
    b.resize(N);
    eigValue.resize(N);
    eigVector.resize(N, N);
    
    
    eigA <<
     1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,10.0,
     2.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,
     3.0,300.0,23.0,24.0,25.0,26.0,27.0,28.0,29.0,30.0,
     4.0,14.0,24.0,34.0,35.0,36.0,37.0,38.0,39.0,40.0,
     5.0,15.0,25.0,35.0,45.0,46.0,47.0,48.0,49.0,50.0,
     6.0,16.0,26.0,36.0,46.0,56.0,57.0,58.0,59.0,60.0,
     7.0,17.0,27.0,37.0,47.0,57.0,67.0,68.0,69.0,70.0,
     8.0,18.0,28.0,38.0,48.0,58.0,68.0,78.0,79.0,80.0,
     9.0,19.0,29.0,39.0,49.0,59.0,69.0,79.0,89.0,90.0,
    10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0;
    
    random_device md;
    
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            A(i, j) = md() % 100;
            eigA(i, j) /= 100.0;
        }
        b(i) = md() % 100;
    }
    
    cout << "行列A(行列式，クラメル用)" << endl;
    cout << A <<  endl;
    cout << "解ベクトル" << endl;
    cout << b <<  endl;
    cout << "行列A(固有値，固有ベクトル用)" << endl;
    cout << eigA << endl;
    //行列式
    tik = clock();
    cout << "行列式(Mydet)：" << MyDet(A) << endl;
    tok = clock();
    cout << "実行時間：" << (double)(tok - tik)/CLOCKS_PER_SEC << "[s]" << endl;
    
    cout << "--------------------------------" <<  endl;
    
    tik = clock();
    cout << "行列式(Eigen)：" << A.determinant() << endl;
    tok = clock();
    cout << "実行時間：" << (double)(tok - tik)/CLOCKS_PER_SEC << "[s]" << endl;
    
    cout << "--------------------------------" <<  endl;
    //
    tik = clock();
    cout << "方程式の解(MyCramer) x = " << endl << MyCramer(A, b) << endl;
    tok = clock();
    cout << "実行時間：" << (double)(tok - tik)/CLOCKS_PER_SEC << "[s]" << endl;
    
    cout << "--------------------------------" <<  endl;
    
    tik = clock();
    cout << "方程式の解(Eigen) x = " << endl << A.fullPivLu().solve(b) << endl;
    tok = clock();
    cout << "実行時間：" << (double)(tok - tik)/CLOCKS_PER_SEC << "[s]" << endl;
    
    cout << "--------------------------------" <<  endl;
    //
    tik = clock();
    eigValue = CalcEigenValue(eigA);
    eigVector = CalcEigenVector(eigA, eigValue);
    cout << "固有値(MyFunc)" << endl << eigValue << endl;
    cout << "固有ベクトル(Myfunc)" << endl << eigVector << endl;
    tok = clock();
    cout << "実行時間：" << (double)(tok - tik)/CLOCKS_PER_SEC << "[s]"  << endl;
    
    cout << "--------------------------------" <<  endl;
    
    tik = clock();
    EigenSolver<Matrix<double,N,N> > es(eigA);
    //SelfAdjointEigenSolver<Matrix<double,N,N>> es(A);
    
    cout << "固有値(Eigen)" << endl << es.eigenvalues() << endl;
    cout << "固有ベクトル(Eigen)" << endl << es.eigenvectors() << endl;
    tok = clock();
    cout << "実行時間：" << (double)(tok - tik)/CLOCKS_PER_SEC << "[s]" << endl;

}

void kadai4_100(){
    clock_t tik, tok;
    MatrixXd A(100,100);
    VectorXd b(100);
    VectorXd eigValue(100);
    MatrixXd eigVector(100, 100);
   
    random_device md;
    
    for(int i = 0; i < 100; i++){
        for(int j = 0; j < 100; j++){
            A(i, j) = md() % 100;
        }
        b(i) = md() % 100;
    }
    
    //行列式
    tik = clock();
    MyDet(A);
    tok = clock();
    cout << "行列式(MyDet)実行時間：" << (double)(tok - tik)/CLOCKS_PER_SEC << "[s]" << endl;
    
    tik = clock();
    A.determinant();
    tok = clock();
    cout << "行列式(Eigen)実行時間：" << (double)(tok - tik)/CLOCKS_PER_SEC << "[s]" << endl;
    
    cout << "--------------------------------" <<  endl;
    //
    tik = clock();
    MyCramer(A, b);
    tok = clock();
    cout << "方程式の解(MyCramer)実行時間：" << (double)(tok - tik)/CLOCKS_PER_SEC << "[s]" << endl;
    
    tik = clock();
    A.fullPivLu().solve(b);
    tok = clock();
    cout << "方程式の解(Eigen)実行時間：" << (double)(tok - tik)/CLOCKS_PER_SEC << "[s]" << endl;
    
    cout << "--------------------------------" <<  endl;
    //
    tik = clock();
    CalcEigenValue(A);
    CalcEigenVector(A, eigValue);
    tok = clock();
    cout << "固有値・固有ベクトル(MyFunc)実行時間：" << (double)(tok - tik)/CLOCKS_PER_SEC << "[s]" << endl;
    
    cout << "--------------------------------" <<  endl;
    
    tik = clock();
    EigenSolver<Matrix<double,100,100> > es(A);
    //es.eigenvalues();
    //es.eigenvectors();
    tok = clock();
    cout << "固有値・固有ベクトル(Eigen)実行時間：" << (double)(tok - tik)/CLOCKS_PER_SEC << "[s]" << endl;
}

int main() {

    kadai1();
    cout << endl;
    kadai2();
    cout << endl;
    kadai3();
    cout << endl;
    kadai4();
    cout << endl;
    //kadai4_100();
    return 0;
}


