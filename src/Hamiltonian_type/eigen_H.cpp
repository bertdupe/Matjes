#include <iostream>
#include <Eigen/Sparse>
#ifdef CPP_MPI
#include <mpi.h>
#endif

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

void eigen_H_bcast(
    int id,
    int mas,
    bool ismas,
    SpMat **Ham,
    MPI_Fint* comm_in){

    MPI::Intracomm comm = MPI_Comm_f2c(*comm_in);
    long shape[3], rows, cols, nnz;
    if(ismas){
        (**Ham).makeCompressed();
        rows=(**Ham).rows();
        cols=(**Ham).cols();
        nnz=(**Ham).nonZeros();
        shape[0]=rows, shape[1]=cols, shape[2]=nnz;
    }
    MPI_Bcast(shape, 3, MPI_LONG,mas, comm);
    rows=shape[0], cols=shape[1], nnz=shape[2];
    if(not ismas){
        *Ham =new SpMat{rows,cols};
        (*Ham)->reserve(nnz);
    }
    int ierr;
    ierr=MPI_Bcast((**Ham).valuePtr(),nnz ,MPI_DOUBLE,mas, comm);
    ierr=MPI_Bcast((**Ham).innerIndexPtr(),nnz ,MPI_INT,mas, comm);
    ierr=MPI_Bcast((**Ham).outerIndexPtr(),cols ,MPI_INT,mas, comm);
    (**Ham).outerIndexPtr()[cols] = nnz;
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

void eigen_H_mult_mat_vec_single(
    SpMat **mat,
    int bnd_min,
    int bnd_max,
    double vec_in[],
    double vec_out[]){

    int cols =(*mat)->cols();
    int size_out = bnd_max-bnd_min+1;
    Map<VectorXd> vec(vec_in,cols);
    Map<VectorXd> res(vec_out,size_out);

    res = ((*mat)->block(bnd_min,0,size_out,cols) * vec);
}


void eigen_H_mult_vec_mat_single(
    SpMat **mat,
    int bnd_min,
    int bnd_max,
    double vec_in[],
    double vec_out[]){

    int rows =(*mat)->rows();
    int size_out = bnd_max-bnd_min+1;
    Map<VectorXd> vec(vec_in,rows);
    Map<VectorXd> res(vec_out,size_out);
    // might be worth checking out if it is better to input vec as colvec instead

    auto slice=(*mat)->block(0,bnd_min,rows,size_out);
    res=slice.transpose() * vec ;
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

void eigen_H_eval_single(
    int ind,
    int dim_mode,
    double vec_l[],
    double vec_r[],
    SpMat **mat,
    double *E){

    int rows=(*mat)->rows();
    Map<Matrix<double,1,Dynamic>> l_vec(vec_l,rows);
    Eigen::SparseVector<double> r_vec(dim_mode);
    for(int i=0; i < dim_mode; i++ ){
        r_vec.coeffRef(i)=vec_r[i];   
    }
    *E = l_vec * (((*mat)->block(0,ind,rows,dim_mode) * r_vec).pruned());
}

void eigen_H_destroy(
    SpMat **Ham){

    delete *Ham;
    *Ham=NULL;
}
}
