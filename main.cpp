#include "MatrixObj.h"
int main(){
    int max  = omp_get_max_threads();
    for (int i = 1; i <= max; ++i) {
        omp_set_num_threads(i);
        std::unique_ptr<MatrixObj> obj(new MatrixObj{512,16, i});
        /*MatrixObj *obj = new MatrixObj{512,16};*/
        obj->setAllMatrices();
        obj->Parallel_runBlockLU_AllInA();
        obj->matrixMultiplication_AllInA();
        obj->countError();
    }




}