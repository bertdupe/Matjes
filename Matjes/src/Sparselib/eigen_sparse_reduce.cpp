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

void reduce_A(const int n,
	      const int m_in,
	      const int m_out,
 	      MatrixXd& A_in,
 	      VectorXd& x_in,
 	      MatrixXd& A_out);


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
