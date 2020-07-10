//compile (decomment main): g++ -std=c++0x  eigen_sparse_reduce.cpp

#include <iostream>
#include <fstream>
#include <iomanip> 
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <math.h>
#include <vector>
#include <Eigen/Core>				
#include <Eigen/SparseCore>
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;


typedef SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of doubles
typedef Triplet<double> Trplt; // declares triplet type of doubles

extern "C"{
void reduce_A(const int n,
	      const int m_in,
	      const int m_out,
 	      MatrixXd& A_in,
 	      VectorXd& x_in,
 	      MatrixXd& A_out);

void get_sparse_matrix(
    double* arr_in,
    const int m_in,
    const int n_in);

//_____________________________________________________________________________________|
//										       |
// 					main (for dev)				       |
//										       |				      
//_____________________________________________________________________________________|
/*
int main(){

	const int n=3;
	const int p=3;
	const int m_in = pow(n,p);  // number of columns
  	const int m_out = pow(n,p-1); 

	// --- declare Eigen types--- //
	VectorXd S_int(n);
	MatrixXd matrix_int(n,m_in);
	MatrixXd matrix_out(n,m_out);

	//---  fill them --- //
	for(int i=0;i<n;i++){
		for(int j=0;j<m_in;j++){
			if(i==j){
				matrix_int(i,j)=i+j;
			}
			else{
				matrix_int(i,j)=0.0;
			}
		}
	}	
	for(int i=0;i<n;i++){
		S_int(i)=i*i;
	}

	cout << "matrix_int=" << endl << matrix_int << endl;
	cout << "S_int=" << endl << S_int << endl;

	// --- call reduce --- //
	 reduce_A(n,m_in,m_out,matrix_int,S_int,matrix_out);

	cout << "matrix_out=" << endl << matrix_out << endl;
}
*/
//_____________________________________________________________________________________|
//										       |
// 					reduce A with eigen			       |
//										       |				      
//_____________________________________________________________________________________|


 //function that uses Eigen to multiply sparse matrix A of size (n,n^p) by dense vector x of size n
void reduce_A(const int n,
	      const int m_in,
	      const int m_out,
 	      MatrixXd& A_in,
 	      VectorXd& x_in,
 	      MatrixXd& A_out){
  
	vector<Trplt> tripletList; // std vector for triplet list: row index, column index, value.
	tripletList.reserve(10*m_in); //reserve room for non-zero elements (here 10/column)

	// --- fill triplet list and fill sparse matrix A--- //
	for(int i=0;i<n;i++){
		for(int j=0;j<m_in;j++){ 
			if (A_in(i,j)!=0.0){ 
	  			tripletList.push_back(Trplt(i,j,A_in(i,j)));
			}
		}	
	}
	SpMat A(n,m_in); //declare sparse matrix of doubles, column-major by default
	A.setFromTriplets(tripletList.begin(), tripletList.end()); //use triplet to fill it
	cout <<" number of non zero elements in A : " << A.nonZeros() << endl; 

	
	// --- get reduced matrix A*x ---//
	VectorXd A_reduced(m_out*n,1); //store in dense matrix 
	A_reduced=A_in.transpose() * x_in;

 	cout <<"reduced matrix has " <<  A_reduced.rows() << " row(s) and " << A_reduced.cols() << " column(s)." << endl;

	// --- reshape into an n,n^{p-1}, ! might need to transpose as well ---//
	Map<MatrixXd> temp(A_reduced.data(), n,m_out);
	A_out=temp;

}


//SparseMatrix<double> all_E;
MatrixXd all_E;


void eigen_set_all_E(
    int m_in,
    int n_in,
    double arr_in[]){
    Map<MatrixXd,RowMajor> m(arr_in,m_in,n_in);
    all_E = m;
    //all_E = m.sparseView();
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
    //VectorXd A= H * vec ;
    //*result = vec.dot(A);
    *result = vec.dot(H * vec);
}

void eigen_matmul_allE(int size_1,double vec_1[],int size_2,double vec_2[],double* result)
{
    Map<RowVectorXd> v1(vec_1,size_1);
    Map<VectorXd> v2(vec_2,size_2);
    *result = v1 * (all_E *v2);
//    VectorXd A= all_E * v2 ;
    //result = v1* (all_E *v2);
}


//void get_sparse_matrix(
//    double arr_in[],
//    const int m_in,
//    const int n_in){
//
//    float data[] = {1,2,3,4};
//    Map<MatrixXf> m(data,2,2);
//
//
// 	cout <<"reduced matrix has " <<  m.rows() << " row(s) and " << m.cols() << " column(s)." << endl;
//
//}

}
