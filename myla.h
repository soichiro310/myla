//
//  myla.h
//  myla
//
//  Created by soichiro310 on 2018/07/21.
//  Copyright © 2018年 soichiro310. All rights reserved.
//

#ifndef myla_h
#define myla_h

#include <iostream>
/*
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/QR>

#include <Eigen/Eigenvalues>
 */
#include <Eigen/Dense>
#include <vector>
#include <cmath>

using namespace std;
using namespace Eigen;

//余因子展開的なノリで行列式を求める
double MyDet(MatrixXd A){
    double det = 1.0;
    double row = A.rows();
    double col = A.cols();
    
    if(row != col || row < 1){
        cerr << "入力された行列は正方行列ではありません" << endl;
        exit(EXIT_FAILURE);
    }
    
    for(int i = 0; i < row; i++){
        for(int j = 0; j < row; j++){
            if (i < j){
                double pivot = A(j, i) / A(i, i);
                for(int k = 0; k < row; k++)
                    A(j, k) -= A(i, k) * pivot;
            }
        }
    }
    
    for(int i = 0; i < row; i++)
        det *= A(i,i);
    
    //cout << "行列式：" << det << endl;
    
    return det;
}

//クラメルの公式による連立方程式
VectorXd MyCramer(MatrixXd A, VectorXd b){
    double detA = MyDet(A);
    double row = A.rows();
    double col = A.cols();
    int size = (int)row;
    VectorXd X;
    X.resize(size);
    
    if(detA == 0){
        cerr << "行列式がゼロとなるため，解が求められません" << endl;
        exit(EXIT_FAILURE);
    }
    
    if(b.rows() != row){
        cerr << "係数行列の行数と，解ベクトルの大きさが同一ではありません" << endl;
        exit(EXIT_FAILURE);
    }
    
    for(int i = 0; i < row; i++){
        MatrixXd Ai = A;
        for(int j = 0; j < col; j++){
            Ai(j, i) = b(j);
        }
        X(i) = MyDet(Ai) / detA;
    }
    /*
    cout << "クラメルの公式による連立方程式の解" << endl;
    for(int i = 0; i < X.size(); i++)
        cout << X(i) << endl;
    */
    return X;
}

//QR法による固有値計算
VectorXd CalcEigenValue(MatrixXd A){
    
    MatrixXd A1;
    A1.resize(A.rows(), A.cols());
    A1 = A;
    for(int i = 0; i < 100; i++){
        HouseholderQR<MatrixXd> qr(A1);
        //ColPivHouseholderQR<MatrixXd> qr(A1);
        MatrixXd Q;
        Q.resize(A.rows(), A.cols());
        Q = qr.householderQ();
        A1 = Q.transpose() * A1 * Q;
    }
    
    VectorXd eigValue;
    int size = (int)A.cols();
    eigValue.resize(size);
    
    for(int i = 0; i < size; i++){
        eigValue(i) = A1(i,i);
    }
    /*
    cout << "QR法による固有値" << endl;
    cout << eigValue << endl;
    */
    return eigValue;
}

//行列とその固有値から，固有ベクトルを求める
MatrixXd CalcEigenVector(MatrixXd A, VectorXd eigValue){
    
    double row = A.rows();
    double col = A.cols();
    MatrixXd B, I, eigV;
    B.resize(row, col);
    I.resize(row, col);
    eigV.resize(row, col);
    I = MatrixXd::Identity(row, col);
    eigV = MatrixXd::Ones(row, col);
    
    for(int i = 0; i < (int)eigValue.size(); i++){
        VectorXd eigVector;
        eigVector.resize(eigValue.size());
        eigVector = VectorXd::Ones(eigValue.size());
        
        B = A - eigValue(i) * I;
        
        for(int j = 0; j < 100; j++){
            eigVector = B.inverse() * eigVector;
            eigVector.normalize();
        }
        
        eigV.col(i) = eigVector;
    }
    /*
    cout << "固有ベクトル" << endl;
    cout << eigV << endl;
    */
    return eigV;
}


#endif /* myla_h */
