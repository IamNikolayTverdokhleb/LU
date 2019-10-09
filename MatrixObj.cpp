#include "MatrixObj.h"
/*Реализация непосредственно алгоритма*/

void MatrixObj::runLU(){
    double t1 = omp_get_wtime();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                U[i] = A[i];
                L[i * n] = A[i * n] / U[0];
                double sum = 0.0;
                    for (int k = 0; k < i; k++){
                         sum += L[i*n+k] * U[k*n+ j];
                    }
                U[i*n + j] = A[i*n+ j] - sum;
                sum = 0.0;
                    if (i > j){
                        L[j*n+i] = 0;
                    }
                         else if(i == j){
                                 L[j*n+i] = 1;
                         }
                             else{
                                for (int k = 0; k < i; k++){
                                          sum += L[j*n+ k] * U[k*n+ i];
                                }
                                    L[j*n + i] = (A[j*n+ i] - sum) / U[i*n + i];
                                }
            }
        }
    double t2 = omp_get_wtime();
    std::cout << "Последовательный алгоритм занял по времени: " << (t2 - t1) << std::endl;
};

void MatrixObj::runBlockLU(){
    double t1 = omp_get_wtime();
/*STEP1: A11 = L11*U11 - basic LU*/
        for (int i = 0; i < r; i++) {
                U[i] = A[i];
                L[i * n] = A[i * n] / U[0];
            for (int j = 0; j < r; j++) {
                    double sum = 0.0;
                    for (int k = 0; k < i; k++){
                        sum += L[i*n+k] * U[k*n+ j];
                    }
                    U[i*n + j] = A[i*n + j] - sum;
                    sum = 0.0;
                    if (i > j){
                         L[j*n+i] = 0;
                         U[n*i+j] = 0.0;
                    }
                        else if(i == j){
                            L[j*n+i] = 1;
                        }
                             else{
                                 for (int k = 0; k < i; k++){
                                     sum += L[j*n+ k] * U[k*n+ i];
                                 }
                                L[j*n + i] = (A[j*n+ i] - sum) / U[i*n + i];
                             }
             }
        }
/*STEP2: A12 = L11*U12 - getting U12*/
    for (int i = 0; i < r; ++i) {
        for (int j = r; j < n; ++j) {
            double sum = 0.0;
                for (int k = 0; k < j; ++k) {
                 sum +=  L[i*n+k]*U[k*n+j];
                 }
            U[i*n+j] = A[i*n+j] - sum;
        }
    }
    /*STEP3: A21 = L21*U11 - getting L21*/
    for (int i = r; i < n; ++i) {
        for (int j = 0; j < r; ++j) {
            double sum = 0;
                for (int k = 0; k < j; k++){
                  sum += L[i*n + k] * U[k*n + j];
                }
            L[i*n + j] = (A[i*n+ j] - sum) / U[j*n + j];
        }
    }
    /*STEP4: A_special22 = A22 - L21*U21 = L22U22 - basic LU*/
    for (int i = r; i < n; ++i) {
        for (int j = r; j < n; ++j) {
            double temp = 0.0;
                for (int k = 0; k < r; ++k) {
                    temp+=L[i*n+k]*U[k*n+j];
                }
            A[i*n+j] -= temp;
        }
    }
    for (int i = 0; i < n-r; i++) {
        for (int j = 0; j < n-r; j++) {
                U[n*r+r+i] = A[n*r+r+i];
                L[n*r+r+i*n] = A[n*r+r+i*n] / U[n*r+r];
            double sum = 0.0;
               for (int k = 0; k < i; k++){
                    sum += L[n*r+r+i*n+k] * U[n*r+r+k*n+ j];
                }
            U[n*r+r+n*i+j] = A[n*r+r+i*n+j] - sum;
            sum = 0.0;
                if (i > j){
                    L[n*r+r+j*n+i] = 0.0;
                    U[n*r+r+n*i+j] = 0.0;
                }
                        else if(i == j){
                            L[n*r+r+j*n+i] = 1;
                        }
                            else{
                                sum = 0;
                                for (int k = 0; k < i; k++){
                                    sum += L[n*r+r+j*n+ k] * U[n*r+r+k*n+ i];
                                }
                                L[n*r+r+j*n + i] = (A[n*r+r+j*n+i] - sum) / U[n*r+r+i*n + i];
                            }
        }
    }
    double t2 = omp_get_wtime();
    std::cout << "Параллельный алгоритм занял по времени: " << (t2 - t1) << std::endl;
};
/*ВОТ ЭТО НЕ РАБОТАЕТ, надо править*/
void MatrixObj::runBlockLUAllInA(){
    double t1 = omp_get_wtime();
/*STEP1: A11 = L11*U11 - basic LU*/
    for (int i = 1; i < r; i++) {
        A[i] = A[i];
        A[i * n] = A[i * n] / A[0];
    }
    for (int i = 1; i < r; i++) {
        for (int j = i; j < r; j++) {
            double sum = 0.0;
                for (int k = 0; k < i; k++){
                    sum += A[i*n+k] * A[k*n+ j];
                }
                A[i*n + j] = A[i*n + j] - sum;
                sum = 0.0;
                for (int k = 0; k < i; k++){
                    sum += A[j*n+ k] * A[k*n+ i];
                }
                A[j*n + i] = (A[j*n+ i] - sum) / A[i*n + i];
        }
    }
    std::cout<<"A MATRIX = " << std::endl;
    getMatrix(A);
/*STEP2: A12 = L11*U12 - getting U12*/
    for (int i = 0; i < r; ++i) {
        for (int j = r; j < n; ++j) {
            double sum = 0.0;
            for (int k = 0; k < j; ++k) {
                sum +=  A[i*n+k]*A[k*n+j];
            }
            A[i*n+j] = A[i*n+j] - sum;
        }
    }
    /*STEP3: A21 = L21*U11 - getting L21*/
    for (int i = r; i < n; ++i) {
        for (int j = 0; j < r; ++j) {
            double sum = 0;
            for (int k = 0; k < j; k++){
                sum += A[i*n + k] * A[k*n + j];
            }
            A[i*n + j] = (A[i*n+ j] - sum) / A[j*n + j];
        }
    }
    /*STEP4: A_special22 = A22 - L21*U21 = L22U22 - basic LU*/
    for (int i = r; i < n; ++i) {
        for (int j = r; j < n; ++j) {
            double temp = 0.0;
            for (int k = 0; k < r; ++k) {
                temp+=A[i*n+k]*A[k*n+j];
            }
            A[i*n+j] -= temp;
        }
    }


    for (int i = 0; i < n-r; i++) {
        A[n*r+r+i] = A[n*r+r+i];
        A[n*r+r+i * n] = A[n*r+r + i * n] / A[0];
    }
    for (int i = 1; i < n-r; i++) {
        for (int j = i; j < n-r; j++) {
            double sum = 0.0;
            for (int k = 0; k < i; k++){
                sum += A[n*r+r+i*n+k] * A[n*r+r+k*n+ j];
            }
            A[n*r+r+i*n + j] = A[n*r+r+i*n + j] - sum;
            sum = 0.0;
            for (int k = 0; k < i; k++){
                sum += A[n*r+r+j*n+ k] * A[n*r+r+k*n+ i];
            }
            A[n*r+r+j*n + i] = (A[n*r+r+j*n+ i] - sum) / A[n*r+r+i*n + i];
        }
    }
    std::cout<<"A MATRIX = " << std::endl;
    getMatrix(A);
    double t2 = omp_get_wtime();
    std::cout << "Параллельный алгоритм занял по времени: " << (t2 - t1) << std::endl;
};
