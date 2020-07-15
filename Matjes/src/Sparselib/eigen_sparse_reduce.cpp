//compile (decomment main): g++ -std=c++0x  -O3 eigen_sparse_reduce.cpp

#include <iostream>
#include <fstream>
#include <iomanip> 
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <math.h>
#include <vector>
#include <time.h>       /* time */
#include <chrono>
#include <Eigen/Core>				
#include <Eigen/SparseCore>
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;


typedef SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of doubles
typedef SparseVector<double> SpVec; // declares a column-major sparse vector type of doubles
typedef Triplet<double> Trplt; // declares triplet type of doubles
typedef vector<double> oneline; //row vector of doubles
typedef vector<oneline> vectoroflines; // vector of rows of doubles


clock_t clock_start, clock_end;
double exe_time;						//execution time											
//_____________________________________________________________________________________|
//										       |
// 					prototypes				       |
//										       |				      
//_____________________________________________________________________________________|

extern "C"{
void build_sparse_matrix(vectoroflines &A_dense,
			 SpMat &A_sparse);

void build_sparse_vector(oneline &x_dense,
			 SpVec &x_sparse);

void reduce_A(const int n,
	      int m_out,
	      SpVec &x_sparse,
	      SpMat &A_sparse_in,
	      SpMat &A_sparse_out);

void reduce_all(const int n,
		const int p,
		SpMat   &A_sparse_in,
		oneline &x_in,
		oneline &A_reduced);


void get_sparse_matrix(
    double* arr_in,
    const int m_in,
    const int n_in);

//_____________________________________________________________________________________|
//										       |
// 					main (for dev)				       |
//										       |				      
//_____________________________________________________________________________________|

int main(){

	const int n=100; //the hamiltonian has dimension n x n^p, this should be passed from Fortran
	const int p=2;

	// --- Eigen types--- //
	 SpMat matrix_sparse_in(n,pow(n,p));

	// ---- std vector types -----//

	vectoroflines matrix_large; 				//hamiltonian
	oneline row;
	oneline M_int;						//order parameter
	oneline matrix_reduced; 				//reduced hamiltonian

	// ---- fill test vectors -----//
	cout << "Filling test matrix of dimensions " << n << " x " << pow(n,p) << " and vector of dim " << n << "..." << endl;
	clock_start = clock();
	for(int i=0;i<n;i++){
		for(int j=0;j<pow(n,p);j++){
			if(i==j){

				row.push_back(i+j);
			}
			else
				row.push_back(0.0);
			}
		matrix_large.push_back(row);
		row.clear();
	}
	
	for(int i=0;i<n;i++){	
		M_int.push_back(i*i);
	}

	/*cout << "matrix_large =" << endl;
	for(int i=0;i<n;i++){
		for(int j=0;j<pow(n,p);j++){
			cout << "\t" << matrix_large[i][j];
		}
		cout << endl;		
	}

	cout << "M_int = " << endl;
	for(int i=0;i<n;i++){
		cout << M_int[i] << endl;
	}*/


	clock_end = clock();
	exe_time = (clock_end - clock_start);
	cout << "Done, execution time = " << exe_time  << " CPS. " << endl;

	// --- build sparse matrix --- //
	cout << "Building sparse matrix..." << endl;
	clock_start = clock();
	build_sparse_matrix(matrix_large,matrix_sparse_in);
	clock_end = clock();
	exe_time = (clock_end - clock_start);
	cout << "Done, execution time = " << exe_time  << " CPS. " << endl;


	// --- call reduce all --- //
	reduce_all(n,p,matrix_sparse_in,M_int,matrix_reduced);
}


//_____________________________________________________________________________________|
//										       |
// 					Build Sparse matrix			       |
//										       |				      
//_____________________________________________________________________________________|

//takes a dense matrix as std vector of rows and outputs eigen sparse matrix

void build_sparse_matrix(vectoroflines &A_dense,
			 SpMat &A_sparse){

	// ----------------- variables---------------------- //

	vector<Trplt> tripletList; 				// std vector for triplet list: row index, column index, value.
	const int n_1=A_dense.size(); 				//number of rows
	const int n_2=A_dense[0].size(); 			//number of columns
	//tripletList.reserve(10* n_2); 				//reserve room for non-zero elements (here 10/column)

	// --- fill triplet list and fill sparse matrix A--- //

	for(int i=0;i<n_1;i++){
		for(int j=0;j<n_2;j++){ 
			if (A_dense[i][j]!=0.0){ 
	  			tripletList.push_back(Trplt(i,j,A_dense[i][j]));
			}
		}	
	}
	A_dense.clear();
	A_sparse.setFromTriplets(tripletList.begin(), tripletList.end()); //use triplet to fill it
 	cout <<"Sparce matrix has "<<  A_sparse.nonZeros() << " non-zero elements."<< endl;
}

//_____________________________________________________________________________________|
//										       |
// 					Build sparse vector			       |
//										       |				      
//_____________________________________________________________________________________|

void build_sparse_vector(oneline &x_dense,
			 SpVec &x_sparse){

	// ----------------- variables---------------------- //


	const int n=x_dense.size(); //number of rows
	x_sparse.reserve(n);
	//vector<Trplt> tripletList; // std vector for triplet list: row index, column index, value

	//tripletList.reserve(n); //reserve room for non-zero elements (here 10/column)

	// --- fill triplet list and fill sparse matrix x--- //

	for(int i=0;i<n;i++){
		if (x_dense[i]!=0.0){ 
	  		x_sparse.insert(i) = x_dense[i];		
		}	
	}
	x_dense.clear();
	//x_sparse.makeCompressed();   
	//x_sparse.setFromTriplets(tripletList.begin(), tripletList.end()); //use triplet to fill it
}
//_____________________________________________________________________________________|
//										       |
// 					reduce A with eigen			       |
//										       |				      
//_____________________________________________________________________________________|

 //function that uses Eigen to multiply sparse matrix A of size (n,n^p) by dense vector x of size n
void reduce_A(const int n,
	      int m_out,
	      SpVec &x_sparse,
	      SpMat &A_sparse_in,
	      SpMat &A_sparse_out){

	// --- perform product ---//
	A_sparse_out=(A_sparse_in.transpose() * x_sparse).pruned();

	// --- reshape into an n,n^{p-1} via dense matrix class---//
	MatrixXd dMat;
	dMat = MatrixXd(A_sparse_out);
	dMat.resize(n,m_out);
	//cout << "dMat="<<endl<<dMat<<endl;
	A_sparse_out=dMat.sparseView();

	//cout << "in reduce_A, after resizing temp=" << endl << MatrixXd(A_sparse_out) << endl;

 	cout <<"reduced matrix has "<<  A_sparse_out.nonZeros() << " non-zero elements, "<<  A_sparse_out.rows() << " row(s) and " << A_sparse_out.cols() << " column(s)." << endl;


}
  
  
  

MatrixXd all_E;

void eigen_set_all_E(
    int m_in,
    int n_in,
    double arr_in[]){
    Map<MatrixXd,RowMajor> m(arr_in,m_in,n_in);
    all_E = m;
    cout <<"set all_E matrix" << endl;
}



SparseMatrix<double> H;
void eigen_set_H(
    int Nentry,
    int Hdim,
    int ind1[],
    int ind2[],
    double arr_in[]){
	typedef Eigen::Triplet<double> T;

	std::vector<T> tripletList;
	tripletList.reserve(Nentry);
	for(int i=0; i < Nentry; i++ ){
		tripletList.push_back(T(ind1[i],ind2[i],arr_in[i]));
	}
	
	SparseMatrix<double> tmp(Hdim,Hdim);
	tmp.setFromTriplets(tripletList.begin(), tripletList.end());
	H=tmp;
 	cout <<"Set H in eigen" << endl;
}
void eigen_eval_H(int dimH,double vec_in[],double* result)
{
    Map<VectorXd> vec(vec_in,dimH);
    *result = vec.dot(H * vec);
}


SparseMatrix<double> B;
void eigen_set_B(
    int Nentry,
    int Hdim,
    int ind1[],
    int ind2[],
    double arr_in[]){
	typedef Eigen::Triplet<double> T;

	std::vector<T> tripletList;
	tripletList.reserve(Nentry);
	for(int i=0; i < Nentry; i++ ){
		tripletList.push_back(T(ind1[i],ind2[i],arr_in[i]));
	}
	
	SparseMatrix<double> tmp(Hdim,Hdim);
	tmp.setFromTriplets(tripletList.begin(), tripletList.end());
	B=tmp;
 	cout <<"Set B in eigen" << endl;
}

void eigen_eval_B(int dimH,double vec_in[],double result[])
{
    Map<VectorXd> vec(vec_in,dimH);
    Map<VectorXd> vec_out(result,dimH);
    vec_out= B * vec ;
}



void eigen_matmul_allE(int size_1,double vec_1[],int size_2,double vec_2[],double* result)
{
    Map<RowVectorXd> v1(vec_1,size_1);
    Map<VectorXd> v2(vec_2,size_2);
    *result = v1 * (all_E *v2);
}

}

//_____________________________________________________________________________________|
//										       |
// 					reduce all				       |
//										       |				      
//_____________________________________________________________________________________|
// takes sparse matrix A and calls reduce_A recursively by doing A=A*x until A is a vector

void reduce_all(const int n,
		const int p, //rank of A -1? so that A is of size n x n^p 
		SpMat   &A_sparse_in,
		oneline &x_in,
		oneline &A_reduced){


	int m_out=pow(n,p-1); 				//number of columns in A after reduction
	SpVec x_sparse(n); 				//sparse vector for order parameter
	SpMat A_sparse_out_1(n,m_out); 			//after one reduction
	SpMat A_sparse_out_2(n,pow(n,p-2)); 		//after two reductions
	MatrixXd A_dense_out;


	//build sparse vector from x
	cout << "Building sparse vector..." << endl;
	clock_start = clock();
	build_sparse_vector(x_in,x_sparse);
	//cout << "x_sparse=" << endl << VectorXd(x_sparse) << endl;
	clock_end = clock();
	exe_time = (clock_end - clock_start);
	cout << "Done, execution time = " << exe_time  << " CPS. " << endl;

	// reduce recursively
	cout << "Beginning reduction..." << endl;
	clock_start = clock();
	reduce_A(n,m_out,x_sparse,A_sparse_in,A_sparse_out_1);
	m_out=pow(n,p-2);
	reduce_A(n,m_out,x_sparse,A_sparse_out_1,A_sparse_out_2);

	//convert back to dense view
	A_dense_out= MatrixXd(A_sparse_out_2);
	clock_end = clock();
	exe_time = (clock_end - clock_start);
	cout << "Done, execution time = " << exe_time  << " CPS. " << endl;
	//cout << "A_dense_out=" << endl << A_dense_out << endl;
	cout << "A_dense_out has " <<A_dense_out.rows() << " row(s) and " <<  A_dense_out.cols() << " column(s). " << endl;

}
