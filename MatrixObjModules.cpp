#include "MatrixObj.h"
/*Реализация всех вспомогательных функций - сеттеров, геттеров, умножений и вычисления погрешности*/
void MatrixObj::setTestMatrix(){
    initialMatrix[0] = 3; initialMatrix[1] = 6; initialMatrix[2] = 8; initialMatrix[3] = 9; initialMatrix[4] = 12;
    initialMatrix[5] = 5; initialMatrix[6] = 11; initialMatrix[7] = 10; initialMatrix[8]= 4; initialMatrix[9] = 8;
    initialMatrix[10] = 7; initialMatrix[11] = 5; initialMatrix[12] = 2; initialMatrix[13] = 1; initialMatrix[14] = 1;
    initialMatrix[15] = 3; initialMatrix[16] = 8; initialMatrix[17] = 0; initialMatrix[18] = -2; A[19] = 4;
    initialMatrix[20] = -5; initialMatrix[21] = 0; initialMatrix[22] = 2; initialMatrix[23] = -1; initialMatrix[24] = 5;
};

void MatrixObj::setAllMatrices(){
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j){
            std::random_device rd;
            std::mt19937 mersenne(rd());
            std::uniform_real_distribution<> urd(0, 1);
            A[i*n + j] = urd(mersenne);;
            initialMatrix[i*n + j] = A[i*n + j];
            U[i*n + j] = 0.0;
            L[i*n + j] = 0.0;
            resMatrix[i*n + j] = 0.0;
        }
};

void MatrixObj::getMatrix(double *Matrix){
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            std::cout << Matrix[i*n + j] << " ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
};

void MatrixObj::getAllMatrices(){
    std::cout<<"A Matrix:"<<std::endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            std::cout << A[i*n + j] << " ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
    std::cout<<"U Matrix:"<<std::endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            std::cout << U[i*n + j] << " ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
    std::cout<<"L Matrix:"<<std::endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            std::cout << L[i*n + j] << " ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
};

void MatrixObj::matrixMultiplication(){
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            for (int k = 0; k < n; ++k) {
                resMatrix[i*n+j] += L[i*n+k]*U[k*n+j];
            }
};

void MatrixObj::matrixMultiplicationAllInA(){
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            for (int k = 0; k < n; ++k) {
                resMatrix[i*n+j] += A[i*n+k]*A[k*n+j];
            }
};

void MatrixObj::countError() {
    double temp = 0.0;
    for (int i = 0; i < n*n; ++i) {
        temp += abs(initialMatrix[i]) - abs(resMatrix[i]);
    }
    std::cout << "Ошибка вычислений равна " << temp << std::endl;
};