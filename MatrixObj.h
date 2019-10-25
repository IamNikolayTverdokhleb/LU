#ifndef LU_MATRIXOBJ_H
#define LU_MATRIXOBJ_H
#include <iostream>
#include <fstream>
#include <iomanip>
#include "omp.h"
#include <random>
class MatrixObj {
public:
    MatrixObj(std::size_t size, std::size_t r, int amount_of_threads):r(r), n(size),amount_of_threads(amount_of_threads), N(size*size){
        std::cout << "Amount of threads = " << amount_of_threads << std::endl;
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
    void setAMatrixFromFile(std::string );
    void getMatrix(double *A) const; /*Геттер для конкретной матрицы*/
    void getAllMatrices() const; /*Геттер для всех матриц*/

    void runLU(); /*Обычное LU разложение*/
    void runLU_AllInA(); /*Обычное LU разложение с записью в матрицу А*/
    void runBlockLU(); /*Блочное LU разложение*/
    void runBlockLU_AllInA(); /*Блочное LU разложение с записью в матрицу А*/
    void Parallel_runBlockLU_AllInA(); /*Параллельное блочное LU разложение с записью в матрицу А*/

    void matrixMultiplication(); /*Перемножение матриц*/
    void matrixMultiplication_AllInA(); /*Перемножение матриц LU, когда они записаны в матрицу A*/
    void countError() const; /*Вычисление ошибки initialMatrix - resMatrix*/
    void splitAMatrixIntoLU(); /*Разбиение матрицы А на L и U*/

private:
    const std::size_t n, N; /*n - размерность строки/столбца матрицы, N - n*n*/
    const std::size_t r; /*размер блока*/
    const int amount_of_threads; /*число потоков для данной реализации*/
    double *A; /*Матрица А*/
    double *L; /*Матрица для L*/
    double *U; /*Матрица для U*/
    double *resMatrix; /*Матрица для результата, получается как перемножение L*U*/
    double *initialMatrix; /*Исходная матрица*/
};


#endif //LU_MATRIXOBJ_H
