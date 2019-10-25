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

void  MatrixObj::runLU_AllInA(){
    double sum = 0.0;
    for (int k = 0; k < n; k++) {
        for (int j = k; j < n; j++) {
                sum = 0;
                for (int p = 0; p < k; p++) {
                    sum += A[k * n + p] * A[p * n + j];
                }
             A[k * n + j] -= sum;
        }

        for (int j = k + 1; j < n; j++) {
                    sum = 0.0;
            for (int p = 0; p < k; p++) {
                     sum += A[j * n + p] * A[p * n + k];
                }
                    A[j * n + k] = (A[j * n + k] - sum) / A[k * n + k];
        }
    }
}

void MatrixObj::runBlockLU_AllInA()
{
    double t1 = omp_get_wtime();
    for (int i = 0; i < n; i += r){
       int  b = (n - i < r) ?  n - i : r;

        /*STEP1: A11 = L11*U11 - basic LU*/
        for (int k = 0; k <= b; k++) {
            for (int j = k + 1; j < b; j++) {
                A[(j + i) * n + k + i] /= A[(k + i) * n + k + i];
                for (int p = k + 1; p < b; p++) {
                    A[(j + i) * n + p + i] -= A[(j + i) * n + k + i] * A[(k + i) * n + p + i];
                }
            }
        }

        int rb = n - i - b; //количество стобцов в U_12 после разложения

        /*STEP2: A12 = L11*U12 - getting U12*/
        for (int j = 1; j < b; j++) {
            for (int k = 0; k < rb; k++) {
                for (int p = 0; p < j; p++) {
                    A[(i + j) * n + i + k + b] -= A[(i + j) * n + i + p] * A[(i + p) * n + i + k + b];
                }
            }
        }
        /*STEP3: A21 = L21*U11 - getting L21*/
        for (int k = 0; k < b; k++) {
            for (int j = 0; j < rb; j++) {
                for (int p = 0; p < k; p++) {
                    A[(i + b + j) * n + i + k] -= A[(i + b + j) * n + i + p] * A[(i + p) * n + i + k];
                }
                A[(i + b + j) * n + i + k] /= A[(i + k) * n + i + k];
            }
        }
        /*STEP4: A_special22 = A22 - L21*U21*/
        for (int p = 0; p < rb; p++)
            for (int k = 0; k < b; k++)
                for (int j = 0; j < rb; j++)
                    A[(p + i + b) * n + (j + i + b)] -= A[(p + i + b) * n + k + i] * A[(k + i) * n + j + i + b];
    }
    double t2 = omp_get_wtime();
    std::cout << "Последовательный блочный алгоритм занял по времени: " << (t2 - t1) << std::endl;
}

void MatrixObj::Parallel_runBlockLU_AllInA(){

    double t1 = omp_get_wtime();
    for (int i = 0; i < n; i += r){
        std::size_t  b = (n - i < r) ?  n - i : r;

        /*STEP1: A11 = L11*U11 - basic LU*/
        for (int k = 0; k <= b; k++) {
            for (int j = k + 1; j < b; j++) {
                A[(j + i) * n + k + i] /= A[(k + i) * n + k + i];
                for (int p = k + 1; p < b; p++) {
                    A[(j + i) * n + p + i] -= A[(j + i) * n + k + i] * A[(k + i) * n + p + i];
                }
            }
        }

        int rb = n - i - b; //количество стобцов в U_12 после разложения

        /*STEP2: A12 = L11*U12 - getting U12*/

        for (int x = 1; x < b; x++) {
            for (int y = 0; y < rb; y++) {
                for (int z = 0; z < x; z++) {
                    A[(i + x) * n + i + y + b] -= A[(i + x) * n + i + z] * A[(i + z) * n + i + y + b];
                }
            }
        }
        /*STEP3: A21 = L21*U11 - getting L21*/
        for (int k = 0; k < b; k++) {
            for (int j = 0; j < rb; j++) {

                for (int p = 0; p < k; p++) {
                    A[(i + b + j) * n + i + k] -= A[(i + b + j) * n + i + p] * A[(i + p) * n + i + k];
                }
                A[(i + b + j) * n + i + k] /= A[(i + k) * n + i + k];
            }
        }
        /*STEP4: A_special22 = A22 - L21*U21*/
        for (int p = 0; p < rb; p++) {
            for (int k = 0; k < b; k++) {
                for (int j = 0; j < rb; j++) {
                    A[(p + i + b) * n + (j + i + b)] -= A[(p + i + b) * n + k + i] * A[(k + i) * n + j + i + b];
                }
            }
        }
    }
    double t2 = omp_get_wtime();
    std::cout << "Параллельный блочный алгоритм занял по времени: " << (t2 - t1) << std::endl;
}