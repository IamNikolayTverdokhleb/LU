#include "MatrixObj.h"
int main(){
    MatrixObj *obj = new MatrixObj{1024,16};

    obj->setAllMatrices();
    obj->Parallel_runBlockLU_AllInA();
    obj->matrixMultiplication_AllInA();
    obj->countError();

}