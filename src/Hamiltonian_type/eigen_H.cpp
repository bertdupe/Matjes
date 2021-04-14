#include <iostream>
#include <Eigen/Sparse>
#include <vector>
#ifdef CPP_MPI
#include <mpi.h>
#endif

using namespace Eigen;

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of doubles
typedef Eigen::SparseVector<double> SpVec; 

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

void eigen_get_transpose(
    SpMat **Ham,
    SpMat **Ham_T){

    int col=(*Ham)->cols();
    int row=(*Ham)->rows();

    *Ham_T =new SpMat{col,row};
  	(**Ham_T)=(*Ham)->transpose().eval();
}


void eigen_H_get_ind_mult_r(
    // find out which indices of the right side (ind_out) are necessary to get all contributions
    // for the ind_in indices of the left vector when applying the Hamiltonian to the right vector
    SpMat **Ham,
    int N_in,       // number of indices of the left vector
    int ind_in[],   // indices of the left vector
    int &N_out,     // in: size of ind_out , out: size of relevant indices in ind_out
    int ind_out[]){

    std::vector<int> out_vec{};
    out_vec.reserve(N_out);
    SpVec vec;
    for (int i=0;i<N_in;i++){
        vec = (**Ham).col(ind_in[i]-1);
        out_vec.insert(out_vec.end(),vec.innerIndexPtr(),vec.innerIndexPtr()+vec.nonZeros());
    }

    std::sort(out_vec.begin(), out_vec.end());
    auto new_end= std::unique(out_vec.begin(), out_vec.end());
    out_vec.erase(new_end, out_vec.end());
    out_vec.shrink_to_fit();

    if(out_vec.size()<(std::vector<int>::size_type) N_out){
        N_out=out_vec.size();
        std::memcpy(ind_out,out_vec.data(),N_out*sizeof(int));
    }
    else{
        std::cerr << "\n\n eigen_H_get_ind_mult_r N_out too small"<< std::endl;
        std::abort();
    }
}

void eigen_H_get_ind_mult_l(
    // find out which indices of the left side (ind_out) are necessary to get all contributions
    // for the ind_in indices of the right vector when applying the Hamiltonian to the left vector
    // slower than eigen_H_get_ind_mult_r
    SpMat **Ham,
    int N_in,       // number of indices of the right vector
    int ind_in[],   // indices of the right vector
    int &N_out,     // in: size of ind_out , out: size of relevant indices in ind_out
    int ind_out[]){

    std::vector<int> out_vec{};
    out_vec.reserve(N_out);
    SpVec vec;
    for (int i=0;i<N_in;i++){
        vec = (**Ham).row(ind_in[i]-1);     // this is slow!
        out_vec.insert(out_vec.end(),vec.innerIndexPtr(),vec.innerIndexPtr()+vec.nonZeros());
    }

    std::sort(out_vec.begin(), out_vec.end());
    auto new_end= std::unique(out_vec.begin(), out_vec.end());
    out_vec.erase(new_end, out_vec.end());
    out_vec.shrink_to_fit();

    if(out_vec.size()<(std::vector<int>::size_type) N_out){
        N_out=out_vec.size();
        std::memcpy(ind_out,out_vec.data(),N_out*sizeof(int));
    }
    else{
        std::cerr << "\n\n eigen_H_get_ind_mult_r N_out too small"<< std::endl;
        std::abort();
    }
}


void eigen_mult_r_disc_disc(
    // slower than eigen_mult_l_disc_disc
    SpMat **Ham,
    int N_r,
    int ind_r[],
    double  vec_r[],
    int N_l,
    int ind_l[],
    double  vec_out[]){

    SpVec r_vec((**Ham).cols());
    r_vec.reserve(N_r); 
    for(int i=0; i < N_r; i++ ){
        r_vec.coeffRef(ind_r[i]-1)=vec_r[i];   
    }
    for(int i=0; i < N_l; i++ ){
        vec_out[i]=(**Ham).row(ind_l[i]-1).dot(r_vec); //this is very slow, as .row() is not contiguous in memory
    }
}

void eigen_mult_l_disc_disc(
    SpMat **Ham,
    int N_l,
    int ind_l[],
    double  vec_l[],
    int N_r,
    int ind_r[],
    double  vec_out[]){

    SpVec l_vec((**Ham).rows());
    l_vec.reserve(N_l); 
    for(int i=0; i < N_l; i++ ){
        l_vec.coeffRef(ind_l[i]-1)=vec_l[i];   
    }
    for(int i=0; i < N_r; i++ ){
        vec_out[i]=(**Ham).col(ind_r[i]-1).dot(l_vec);
    }
}


#ifdef CPP_MPI


void eigen_H_send(
    const int& id,
    const int& tag,
    SpMat **Ham,
    MPI_Fint* comm_in){

    MPI::Intracomm comm = MPI_Comm_f2c(*comm_in);
    long shape[3], rows, cols, nnz;
    (**Ham).makeCompressed();
    rows=(**Ham).rows();
    cols=(**Ham).cols();
    nnz=(**Ham).nonZeros();
    shape[0]=rows, shape[1]=cols, shape[2]=nnz;

    int ierr;
    ierr=MPI_Send(shape, 3, MPI_LONG, id, tag,  comm);
    ierr=MPI_Send((**Ham).valuePtr(),      nnz,  MPI_DOUBLE, id, tag, comm);
    ierr=MPI_Send((**Ham).innerIndexPtr(), nnz,  MPI_INT,    id, tag, comm);
    ierr=MPI_Send((**Ham).outerIndexPtr(), cols, MPI_INT,    id, tag, comm);
}

void eigen_H_recv(
    const int& id,
    const int& tag,
    SpMat **Ham,
    MPI_Fint* comm_in){

    MPI::Intracomm comm = MPI_Comm_f2c(*comm_in);
    MPI_Status status;
    int ierr;
    long shape[3], rows, cols, nnz;

    
    ierr=MPI_Recv(shape, 3, MPI_LONG, id, tag,  comm, &status);
    rows=shape[0], cols=shape[1], nnz=shape[2];
    *Ham =new SpMat{rows,cols};
    (*Ham)->reserve(nnz);
    ierr=MPI_Recv((**Ham).valuePtr(),      nnz,  MPI_DOUBLE, id, tag, comm, &status);
    ierr=MPI_Recv((**Ham).innerIndexPtr(), nnz,  MPI_INT,    id, tag, comm, &status);
    ierr=MPI_Recv((**Ham).outerIndexPtr(), cols, MPI_INT,    id, tag, comm, &status);
    (**Ham).outerIndexPtr()[cols] = nnz;
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
#endif

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

void eigen_H_mult_mat_vec_cont(
    SpMat **mat,
    int bnd_min,
    int bnd_max,
    double vec_in[],
    double vec_out[]){

    int size_in = bnd_max-bnd_min+1;
    Map<VectorXd> res(vec_out,(**mat).rows());

    SpVec vec((**mat).cols());
    vec.reserve(size_in);
    for(int i=0; i < size_in; i++ ){
        vec.coeffRef(bnd_min+i)=vec_in[i];   
    }

    res = **mat * vec;
}

void eigen_H_mult_vec_mat_cont(
    SpMat **mat,
    int bnd_min,
    int bnd_max,
    double vec_in[],
    double vec_out[]){

    int size_in = bnd_max-bnd_min+1;
    Map<VectorXd> res(vec_out,(**mat).cols());

    SpVec vec((**mat).rows());
    vec.reserve(size_in);
    for(int i=0; i < size_in; i++ ){
        vec.coeffRef(bnd_min+i)=vec_in[i];   
    }

    res= (*mat)->transpose() * vec ;
}


void eigen_H_mult_mat_vec_disc(
    SpMat **mat,
    int N,
    int ind[0],
    double val[],
    double vec_out[]){

    Map<VectorXd> res(vec_out,(**mat).rows());
    SpVec vec((**mat).cols());
    vec.reserve(N);
    for(int i=0; i < N ; i++ ){
        vec.coeffRef(ind[i]-1)=val[i];   
    }
    res = **mat * vec;
}

void eigen_H_mult_vec_mat_disc(
    SpMat **mat,
    int N,
    int ind[0],
    double val[],
    double vec_out[]){

    Map<VectorXd> res(vec_out,(**mat).cols());
    SpVec vec((**mat).rows());
    vec.reserve(N);
    for(int i=0; i < N ; i++ ){
        vec.coeffRef(ind[i]-1)=val[i];   
    }
    res= (**mat).transpose() * vec ;
}


void eigen_H_mult_l_ind(
    SpMat **mat,
    double vec_in[],
    int& N,
    int ind_out[],
    double vec_out[]){

    Map<VectorXd> vec(vec_in,(*mat)->rows());
    for(std::vector<SpMat>::size_type i=0; i < (std::vector<SpMat>::size_type) N ; i++ ){
        vec_out[i]=(**mat).col(ind_out[i]-1).dot(vec);
    }
}

void eigen_H_mult_r_ind(
    SpMat **mat,
    double vec_in[],
    int& N,
    int ind_out[],
    double vec_out[]){

    Map<VectorXd> vec(vec_in,(*mat)->cols());
    for(std::vector<SpMat>::size_type i=0; i < (std::vector<SpMat>::size_type) N ; i++ ){
        vec_out[i]=(**mat).row(ind_out[i]-1).dot(vec);
    }
}

void eigen_H_mult_mat_disc_disc(
    SpMat **mat,
    int N_in,
    int ind_in[0],
    double vec_in[],
    int& N_out,  // in: size of ind_out , out: size of relevant indices in ind_out
    int ind_out[],
    double vec_out[]){

    int rows =(*mat)->rows();

    SpVec vec((**mat).cols());
    vec.reserve(N_in);
    for(int i=0; i < N_in ; i++ ){
        vec.coeffRef(ind_in[i]-1)=vec_in[i];
    }
    SpVec res(rows);
    res= (**mat) * vec ;
    if(res.nonZeros()<(std::vector<int>::size_type) N_out){
        N_out=res.nonZeros();
        std::memcpy(vec_out,res.valuePtr(),N_out*sizeof(double));
        std::memcpy(ind_out,res.innerIndexPtr(),N_out*sizeof(int));
        for(int i=0;i < N_out;i++) ind_out[i]++;
    }
    else{
        std::cerr << "\n\n eigen_H_mult_mat_disc_disc N_out too small"<< std::endl;
        std::cerr << "initial N_out "<< N_out<<std::endl;
        std::cerr << "size output   "<< res.nonZeros()<<std::endl;
        std::abort();
    }

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
    SpVec r_vec(dim_mode);
    for(int i=0; i < dim_mode; i++ ){
        r_vec.coeffRef(i)=vec_r[i];   
    }
    *E = l_vec * (((*mat)->block(0,ind,rows,dim_mode) * r_vec).pruned());
}

void eigen_H_destroy(
    SpMat **Ham){

    delete *Ham;
    Ham=NULL;
}
}
