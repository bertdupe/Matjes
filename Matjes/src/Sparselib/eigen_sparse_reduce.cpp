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
typedef vector<double> oneline; //row vector of doubles
typedef vector<oneline> vectoroflines; // vector of rows of doubles


//_____________________________________________________________________________________|
//										       |
// 					prototypes				       |
//										       |				      
//_____________________________________________________________________________________|


void build_sparse_matrix(vectoroflines &A_dense,
			 SpMat &A_sparse);

void reduce_A(const int n,
	      int m_out,
	      oneline &x_in,
	      SpMat &A_sparse_in,
	      SpMat &A_sparse_out);

void reduce_all(const int p, 
		SpMat   &A_sparse_in,
		oneline &x_in,
		oneline &A_reduced);



//_____________________________________________________________________________________|
//										       |
// 					main (for dev)				       |
//										       |				      
//_____________________________________________________________________________________|

int main(){

	const int n=3;
	const int p=3;


	// --- declare Eigen types--- //
	 SpMat matrix_sparse_in;

	// ---- std vector types -----//

	vectoroflines matrix_large; //hamiltonian
	oneline row;
	oneline M_int; //order parameter
	oneline matrix_reduced; //reduced hamiltonian

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

	cout << "matrix_large =" << endl;
	for(int i=0;i<n;i++){
		for(int j=0;j<pow(n,p);j++){
			cout << "\t" << matrix_large[i][j];
		}
		cout << endl;		
	}

	cout << "M_int = " << endl;
	for(int i=0;i<n;i++){
		cout << M_int[i] << endl;
	}


	// --- build sparse matrix --- //
	build_sparse_matrix(matrix_large,matrix_sparse_in);


	// --- call reduce all --- //
	reduce_all(p,matrix_sparse_in,M_int,matrix_reduced);

	//cout << "matrix_out=" << endl << matrix_out << endl;
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

	vector<Trplt> tripletList; // std vector for triplet list: row index, column index, value.
	//tripletList.reserve(10*m_in); //reserve room for non-zero elements (here 10/column)
	const int n_1=A_dense.size(); //number of rows
	const int n_2=A_dense[0].size(); //number of columns

	// --- fill triplet list and fill sparse matrix A--- //

	for(int i=0;i<n_1;i++){
		for(int j=0;j<n_2;j++){ 
			if (A_dense[i][j]!=0.0){ 
	  			tripletList.push_back(Trplt(i,j,A_dense[i][j]));
			}
		}	
	}
	//SpMat A_sparse(n_1,n_2); //declare sparse matrix of doubles, column-major by default
	A_sparse.setFromTriplets(tripletList.begin(), tripletList.end()); //use triplet to fill it
}

//_____________________________________________________________________________________|
//										       |
// 					reduce A with eigen			       |
//										       |				      
//_____________________________________________________________________________________|

 //function that uses Eigen to multiply sparse matrix A of size (n,n^p) by dense vector x of size n
void reduce_A(const int n,
	      int m_out,
	      oneline &x_in,
	      SpMat &A_sparse_in,
	      SpMat &A_sparse_out){

	//map std vector to eigen type
    	VectorXd x =Map< Matrix<double, Dynamic, 1> >(x_in.data());
	//SpMat x= Map<SpMat,Dynamic,1>(x_in.data());
  	// Map<MatrixXd> temp(A_reduced.data(), n,m_out);
	// --- get reduced matrix A*x ---//
	//VectorXd A_reduced(m_out*n,1); //store in dense matrix 
	//SpMat A_reduced(m_out*n,1); //store in sparse matrix 	
	A_sparse_out=(A_sparse_in.transpose() * x).pruned(); 
	//A_sparse_out.prune();

//sm3 = (sm1 * sm2).pruned();                  // product of sparse matrices removes numerical zeros
//mat.prune()

 	cout <<"reduced matrix has " A_sparse_out.nonZeros() << " non-zero elements, "<<  A_sparse_out.rows() << " row(s) and " << A_sparse_out.cols() << " column(s)." << endl;

	// --- reshape into an n,n^{p-1}, ! might need to transpose or not ---//
	//SparseMatrix<double> A11;
   	 A_sparse_out.resize(n,m_out);

	//Map<MatrixXd> temp(A_reduced.data(), n,m_out);
	//A_out=temp.transpose();

}

//_____________________________________________________________________________________|
//										       |
// 					reduce all				       |
//										       |				      
//_____________________________________________________________________________________|
// takes sparse matrix A and calls reduce_A recursively by doing A=A*x until A is a vector

void reduce_all(const int p, //rank of A -1? so that A is of size n x n^p 
		SpMat   &A_sparse_in,
		oneline &x_in,
		oneline &A_reduced){

	const int n=x_in.size(); 
	//const int m_in=A_dense[0].size(); //number of columns in A	
	int m_out=pow(n,p-1); //number of columns in A after reduction

	//SpMat A_sparse_in;
	SpMat A_sparse_out; 
	//MatrixXd A_dense_out(dynamical?); //output of reduce_A, size should be dynamical


	//reduce A_sparse until it is a vector
	//while(m_out>1){
		//reduce_A(n,m_out,A_sparse_in,x_in,A_sparse_out);
	reduce_A(n,m_out,x_in,A_sparse_in,A_sparse_out);

	//	m_out=A_dense_out.cols()

		//then something like A_sparse = A_dense_out 

	//}


}

