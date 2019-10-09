#ifndef LU_MATRIXOBJ_H
#define LU_MATRIXOBJ_H
#include <iostream>
#include <iomanip>
#include "omp.h"
#include <random>
class MatrixObj {
public:
    MatrixObj(std::size_t size, std::size_t r){
        this->r =  r;
        this->n = size;
        this->N = size*size;
        A = new double[N];
        L = new double[N];
        U = new double[N];
        resMatrix = new double[N];
        initialMatrix = new double[N];
    };
    ~MatrixObj(){
        delete[] A;
        delete[] L;
        delete[] U;
        delete[] initialMatrix;
        delete[] resMatrix;
    };
    void setTestMatrix(); /*Задание тестовой матрицы*/
    void setAllMatrices(); /*Задание рандомной матрицы А = resMatrix и нулевых L,U*/
    void getMatrix(double *A); /*Геттер для конкретной матрицы*/
    void getAllMatrices(); /*Геттер для всех матриц*/
    void runLU(); /*Обычное LU разложение*/
    void runBlockLU(); /*Блочное LU разложение*/
    void runBlockLUAllInA(); /*Блочное LU разложение с записью в матрицу А*/
    void matrixMultiplication(); /*Перемножение матриц*/
    void matrixMultiplicationAllInA(); /*Перемножение матриц LU, когда они записаны в матрицу A*/
    void countError(); /*Вычисление ошибки initialMatrix - resMatrix*/
private:
    std::size_t n, N; /*n - размерность строки/столбца матрицы, N - n*n*/
    std::size_t r; /*размер блока*/
    double *A; /*Матрица А*/
    double *L; /*Матрица для L*/
    double *U; /*Матрица для U*/
    double *resMatrix; /*Матрица для результата, получается как перемножение L*U*/
    double *initialMatrix; /*Исходная матрица*/
};


#endif //LU_MATRIXOBJ_H
