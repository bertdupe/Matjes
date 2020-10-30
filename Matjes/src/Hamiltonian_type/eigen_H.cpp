#include <iostream>
#include <Eigen/Sparse>

using namespace Eigen;

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of doubles

extern "C"{

void eigen_H_init(
    int Nentry,
    int Hdim[2],
    int rows[],
    int cols[],
    double arr_in[],
    SpMat **Ham){

	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(Nentry);
	for(int i=0; i < Nentry; i++ ){
		tripletList.push_back(T(rows[i],cols[i],arr_in[i]));
	}

    *Ham =new SpMat{Hdim[0],Hdim[1]};
  	(*Ham)->setFromTriplets(tripletList.begin(), tripletList.end());
}

void eigen_H_add(
    SpMat **Ham_sum,
    SpMat **Ham_add){

    **Ham_sum=**Ham_sum+**Ham_add;
}

void eigen_H_mult_mat_vec(
    SpMat **mat,
    double vec_in[],
    double vec_out[]){

    Map<VectorXd> vec(vec_in,(*mat)->cols());
    Map<VectorXd> res(vec_out,(*mat)->rows());
    res= (**mat) * vec;
}

void eigen_H_mult_vec_mat(
    SpMat **mat,
    double vec_in[],
    double vec_out[]){

    Map<VectorXd> vec(vec_in,(*mat)->rows());
    Map<VectorXd> res(vec_out,(*mat)->cols());
    res=(*mat)->transpose() * vec ;
}

void eigen_H_copy(
    SpMat **Ham_in,
    SpMat **Ham_out){

    *Ham_out=new SpMat{(*Ham_in)->rows(),(*Ham_in)->cols()};
    **Ham_out=**Ham_in;
}

void eigen_H_destroy(
    SpMat **Ham){

    delete *Ham;
    *Ham=NULL;
}
}
