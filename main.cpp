#include "MatrixObj.h"
int main() {
    /*Указываем размер матрицы(n) и блока (r)*/
    std::size_t n = 10;
    std::size_t r = 7;

    MatrixObj *obj = new MatrixObj(n,r); /*Создаем обьект*/
    obj->setAllMatrices(); /*Заполняем все матрицы*/
    obj->runBlockLU(); /*Запускаем блочное LU*/
    obj->matrixMultiplication(); /*LU = resMatrix*/
    obj->countError(); /*initialMatrix - resMatrix*/
}