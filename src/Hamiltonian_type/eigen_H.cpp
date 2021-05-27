#include <iostream>
#include <Eigen/Sparse>
#include <vector>
#include <numeric>

#ifdef CPP_MPI
#include <mpi.h>
#endif

using namespace Eigen;

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of doubles
typedef Eigen::SparseVector<double> SpVec; 

extern "C"{

void eigen_H_init(
    int Nentry,
    const int Hdim[2],
    const int rows[],
    const int cols[],
    const double arr_in[],
    SpMat* &Ham){

	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(Nentry);
	for(int i=0; i < Nentry; i++ ){
		tripletList.push_back(T(rows[i],cols[i],arr_in[i]));
	}

    Ham =new SpMat{Hdim[0],Hdim[1]};
  	Ham->setFromTriplets(tripletList.begin(), tripletList.end());
}

void eigen_get_transpose(
    const SpMat* &Ham,
    SpMat* &Ham_T){

    int col=Ham->cols();
    int row=Ham->rows();

    Ham_T =new SpMat{col,row};
  	*Ham_T=Ham->transpose().eval();
}

#ifdef CPP_MPI
void eigen_H_send(
    const int& id,
    const int& tag,
    const SpMat* &Ham,
    const MPI_Fint* comm_in){

    MPI::Intracomm comm = MPI_Comm_f2c(*comm_in);
    long shape[3], rows, cols, nnz;
    Ham->makeCompressed();
    rows=Ham->rows();
    cols=Ham->cols();
    nnz=Ham->nonZeros();
    shape[0]=rows, shape[1]=cols, shape[2]=nnz;

    int ierr;
    ierr=MPI_Send(shape, 3, MPI_LONG, id, tag,  comm);
    ierr=MPI_Send(Ham->valuePtr(),      nnz,  MPI_DOUBLE, id, tag, comm);
    ierr=MPI_Send(Ham->innerIndexPtr(), nnz,  MPI_INT,    id, tag, comm);
    ierr=MPI_Send(Ham->outerIndexPtr(), cols, MPI_INT,    id, tag, comm);
}

void eigen_H_recv(
    const int& id,
    const int& tag,
    SpMat* &Ham,
    const MPI_Fint* comm_in){

    MPI::Intracomm comm = MPI_Comm_f2c(*comm_in);
    MPI_Status status;
    int ierr;
    long shape[3], rows, cols, nnz;

    
    ierr=MPI_Recv(shape, 3, MPI_LONG, id, tag,  comm, &status);
    rows=shape[0], cols=shape[1], nnz=shape[2];
    Ham =new SpMat{rows,cols};
    Ham->reserve(nnz);
    ierr=MPI_Recv(Ham->valuePtr(),      nnz,  MPI_DOUBLE, id, tag, comm, &status);
    ierr=MPI_Recv(Ham->innerIndexPtr(), nnz,  MPI_INT,    id, tag, comm, &status);
    ierr=MPI_Recv(Ham->outerIndexPtr(), cols, MPI_INT,    id, tag, comm, &status);
    Ham->outerIndexPtr()[cols] = nnz;
}


void eigen_H_bcast(
    const int &id,
    const int &mas,
    const bool &ismas,
    SpMat* &Ham,
    const MPI_Fint* comm_in){

    MPI::Intracomm comm = MPI_Comm_f2c(*comm_in);
    long shape[3], rows, cols, nnz;
    if(ismas){
        Ham->makeCompressed();
        rows=Ham->rows();
        cols=Ham->cols();
        nnz= Ham->nonZeros();
        shape[0]=rows, shape[1]=cols, shape[2]=nnz;
    }
    MPI_Bcast(shape, 3, MPI_LONG,mas, comm);
    rows=shape[0], cols=shape[1], nnz=shape[2];
    if(not ismas){
        Ham =new SpMat{rows,cols};
        Ham->reserve(nnz);
    }
    int ierr;
    ierr=MPI_Bcast(Ham->valuePtr(),nnz ,MPI_DOUBLE,mas, comm);
    ierr=MPI_Bcast(Ham->innerIndexPtr(),nnz ,MPI_INT,mas, comm);
    ierr=MPI_Bcast(Ham->outerIndexPtr(),cols ,MPI_INT,mas, comm);
    Ham->outerIndexPtr()[cols] = nnz;
}

#if 1
void eigen_H_distribute(
    const int& id,
    const int& mas,
    const bool& ismas,
    SpMat* &Ham,
    const MPI_Fint* comm_in){

    MPI::Intracomm comm = MPI_Comm_f2c(*comm_in);
    long shape[3], rows, cols, nnz;
    int ierr;
    if(ismas){
        Ham->makeCompressed();
        rows=Ham->rows();
        cols=Ham->cols();
        nnz=Ham->nonZeros();
        shape[0]=rows, shape[1]=cols, shape[2]=nnz;
    }
    MPI_Bcast(shape, 3, MPI_LONG,mas, comm);
    rows=shape[0], cols=shape[1], nnz=shape[2];
    int Np;
    MPI_Comm_size(comm, &Np);

    //get array with aim number of entries per matrix
    std::vector<int> nnz_count;
    nnz_count.resize(Np);
    std::fill(nnz_count.begin(), nnz_count.end(), nnz/Np);
    for( int i = 0; i<nnz%Np; i++) nnz_count[i]++;
    std::vector<int> displ;
    displ.resize(Np);
    std::partial_sum(nnz_count.begin(), nnz_count.end()-1, displ.begin()+1);

    if(not ismas){
        *Ham =new SpMat{rows,cols};
        (*Ham)->reserve(nnz_count[id]);
    }
    ierr=MPI_Bcast(Ham->outerIndexPtr(),cols ,MPI_INT,mas, comm);
    auto point=Ham->outerIndexPtr();
    for(int i=0;i<cols; i++){
        point[i]=point[i]-displ[id];
        point[i]=std::max(point[i],0);
        point[i]=std::min(point[i],nnz_count[id]);
    }

    if(ismas){
        SpMat* tmp =new SpMat{rows,cols};
        (tmp)->reserve(nnz_count[id]);
        ierr=MPI_Scatterv(Ham->innerIndexPtr(), nnz_count.data(), displ.data(), MPI_INT,    (*tmp).innerIndexPtr(), nnz_count[id], MPI_INT,    mas, comm);
        ierr=MPI_Scatterv(Ham->valuePtr(),      nnz_count.data(), displ.data(), MPI_DOUBLE, (*tmp).valuePtr(),      nnz_count[id], MPI_DOUBLE, mas, comm);
        std::memcpy((*tmp).outerIndexPtr(),Ham->outerIndexPtr(),cols*sizeof(Eigen::SparseMatrix<double>::StorageIndex));
        delete *Ham;
        *Ham=tmp;
    }
    else{
        ierr=MPI_Scatterv(NULL,                    nnz_count.data(), displ.data(), MPI_INT,    Ham->innerIndexPtr(), nnz_count[id], MPI_INT,    mas, comm);
        ierr=MPI_Scatterv(NULL,                    nnz_count.data(), displ.data(), MPI_DOUBLE, Ham->valuePtr(),      nnz_count[id], MPI_DOUBLE, mas, comm);
    }
    Ham->outerIndexPtr()[cols] = nnz_count[id];
    Ham->makeCompressed();
}
#else
void eigen_H_distribute(
    // intermediate csc state which distributes the matrices differently ( fewer entries per row instead of empty rows)
    const int& id,
    const int& mas,
    const bool& ismas,
    SpMat* &Ham,
    MPI_Fint* comm_in){

    MPI::Intracomm comm = MPI_Comm_f2c(*comm_in);
    long shape[3], rows, cols, nnz;
    int ierr;
    if(ismas){
        Ham->makeCompressed();
        rows=Ham->rows();
        cols=Ham->cols();
        nnz=Ham->nonZeros();
        shape[0]=rows, shape[1]=cols, shape[2]=nnz;
    }
    MPI_Bcast(shape, 3, MPI_LONG,mas, comm);
    rows=shape[0], cols=shape[1], nnz=shape[2];
    int Np;
    MPI_Comm_size(comm, &Np);

    //get array with aim number of entries per matrix
    std::vector<int> nnz_count;
    nnz_count.resize(Np);
    std::fill(nnz_count.begin(), nnz_count.end(), nnz/Np);
    for( int i = 0; i<nnz%Np; i++) nnz_count[i]++;
    std::vector<int> displ;
    displ.resize(Np);
    std::partial_sum(nnz_count.begin(), nnz_count.end()-1, displ.begin()+1);

    Eigen::SparseMatrix<double,RowMajor> mat_row{rows,cols}; 
    mat_row.reserve(nnz_count[id]);
    Eigen::SparseMatrix<double,RowMajor> mat_tmp; 
    if(ismas){
        mat_tmp=Eigen::SparseMatrix<double,RowMajor>((*Ham));
        delete *Ham;
        std::memcpy(mat_row.outerIndexPtr(),mat_tmp.outerIndexPtr(),cols*sizeof(Eigen::SparseMatrix<double>::StorageIndex));
    }

    ierr=MPI_Bcast(mat_row.outerIndexPtr(),cols ,MPI_INT,mas, comm);

    auto point=(mat_row).outerIndexPtr();
    for(int i=0;i<cols; i++){
        point[i]=point[i]-displ[id];
        point[i]=std::max(point[i],0);
        point[i]=std::min(point[i],nnz_count[id]);
    }

    if(ismas){
        ierr=MPI_Scatterv(mat_tmp.innerIndexPtr(), nnz_count.data(), displ.data(), MPI_INT,    mat_row.innerIndexPtr(), nnz_count[id], MPI_INT,    mas, comm);
        ierr=MPI_Scatterv(mat_tmp.valuePtr(),      nnz_count.data(), displ.data(), MPI_DOUBLE, mat_row.valuePtr(),      nnz_count[id], MPI_DOUBLE, mas, comm);
    }
    else{
        ierr=MPI_Scatterv(NULL,                    nnz_count.data(), displ.data(), MPI_INT,    mat_row.innerIndexPtr(), nnz_count[id], MPI_INT,    mas, comm);
        ierr=MPI_Scatterv(NULL,                    nnz_count.data(), displ.data(), MPI_DOUBLE, mat_row.valuePtr(),      nnz_count[id], MPI_DOUBLE, mas, comm);
    }
    mat_row.outerIndexPtr()[rows] = nnz_count[id];
    mat_row.makeCompressed();

    *Ham =new SpMat(mat_row);
    Ham->makeCompressed();
}
#endif

#endif

void eigen_H_add(
    SpMat* &Ham_sum,
    const SpMat* &Ham_add){

    *Ham_sum=*Ham_sum+*Ham_add;
}

void eigen_H_mult_r(
    const SpMat* &mat,
    const double vec_in[],
    double vec_out[]){

    Map<const VectorXd> vec(vec_in,mat->cols());
    Map<VectorXd> res(vec_out,mat->rows());
    res= (*mat) * vec;
}

void eigen_H_mult_r_beta(
    const SpMat*  &mat,
    const double  vec_in[],
    double        vec_out[],
    const double  &beta){

    Map<const VectorXd> vec(vec_in,mat->cols());
    Map<VectorXd> res(vec_out,mat->rows());
    res= (*mat) * vec + beta * res;
}


void eigen_H_mult_l_ind(
    const SpMat* &mat,
    const double vec_in[],
    const int& N,
    const int ind_out[],
    double vec_out[]){

    Map<const VectorXd> vec(vec_in,mat->rows());
    for(std::vector<SpMat>::size_type i=0; i < (std::vector<SpMat>::size_type) N ; i++ ){
        vec_out[i]=mat->col(ind_out[i]-1).dot(vec);
    }
}

void eigen_H_mult_r_ind(
    const SpMat* &mat,
    const double vec_in[],
    const int& N,
    const int ind_out[],
    double vec_out[]){

    Map<const VectorXd> vec(vec_in,mat->cols());
    for(std::vector<SpMat>::size_type i=0; i < (std::vector<SpMat>::size_type) N ; i++ ){
        vec_out[i]=mat->row(ind_out[i]-1).dot(vec);
    }
}

void eigen_get_dat(
    //pass the data array pointers to the outside
    SpMat* &mat,
    int &nnz,
    int size[2],
    int* &col,
    int* &row,
    double* &val){
    if(! mat->isCompressed()) mat->makeCompressed();
    nnz=mat->nonZeros();
    size[0]=mat->rows();
    size[1]=mat->cols();
    val=mat->valuePtr();
    row=mat->innerIndexPtr();
    col=mat->outerIndexPtr();
}


void eigen_H_mult_l(
    const SpMat* &mat,
    const double vec_in[],
    double vec_out[]){

    Map<const VectorXd> vec(vec_in,mat->rows());
    Map<VectorXd> res(vec_out,mat->cols());
    res=mat->transpose() * vec;
}

void eigen_H_mult_l_beta(
    const SpMat* &mat,
    const double vec_in[],
    double vec_out[],
    const double  &beta){

    Map<const VectorXd> vec(vec_in,mat->rows());
    Map<VectorXd> res(vec_out,mat->cols());
    res=mat->transpose() * vec + beta * res;
}

void eigen_H_copy(
    const SpMat* &Ham_in,
    SpMat* &Ham_out){

    Ham_out=new SpMat{Ham_in->rows(),Ham_in->cols()};
    *Ham_out=*Ham_in;
}

void eigen_H_destroy(
    SpMat* &Ham){

    delete Ham;
    Ham=NULL;
}
}
