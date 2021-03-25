#include <iostream>
#include <Eigen/Sparse>
#include <vector>

using namespace Eigen;

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of doubles

extern "C"{


void modes_alloc(
    int N,
    std::vector<SpMat> **modes){

    *modes =new std::vector<SpMat>;
    (*modes)->resize(N);
}

void modes_init(
    std::vector<SpMat> **modes,
    int i_mode,
    int nnz,
    int dim[2],
    int rows[],
    int cols[],
    double val[]){

	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(nnz);
	for(int i=0; i < nnz; i++ ){
		tripletList.push_back(T(rows[i],cols[i],val[i]));
	}

    SpMat* mode=new SpMat{dim[0],dim[1]};
    mode->setFromTriplets(tripletList.begin(), tripletList.end());
    (*modes)->at(i_mode)=(*mode);
}

void get_mode_i(
    std::vector<SpMat> **modes,
    int i_mode,
    double vec_in[],
    double vec_out[]){
    Map<VectorXd> vec(vec_in ,(*modes)->at(i_mode).cols());
    Map<VectorXd> res(vec_out,(*modes)->at(i_mode).rows());

    res= (*modes)->at(i_mode)* vec ;
}

void modes_set(
    int ind,
    std::vector<SpMat> *modes){

    modes =new std::vector<SpMat>;
    modes->reserve(ind);
}

void modes_copy(
    std::vector<SpMat> **modes_in,
    std::vector<SpMat> **modes_out){

    *modes_out =new std::vector<SpMat>;
    (*modes_out)->resize((*modes_in)->size());
    int rows, cols;
    for(std::vector<SpMat>::size_type i=0; i<(*modes_in)->size();i++) {
        rows=(*modes_in)->at(i).rows();
        cols=(*modes_in)->at(i).cols();
        SpMat* mode=new SpMat{rows,cols};
        *mode=(*modes_in)->at(i);
        (*modes_out)->at(i)=(*mode);
    }
}

void mode_reduce(
    std::vector<SpMat> **modes,
    int i_mode,
    double vec_in[],
    double vec_out[]){

    Map<VectorXd> in(vec_in,(*modes)->at(i_mode).rows());
    Map<VectorXd> out(vec_out,(*modes)->at(i_mode).cols());
    out=(*modes)->at(i_mode).transpose() * in ;
}

void mode_destroy(
    std::vector<SpMat> **modes){
    
    (**modes).clear();
    delete *modes;
    modes=NULL;
}


}
