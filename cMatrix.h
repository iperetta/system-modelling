/* code by: iperetta@ieee.org */
#ifndef MYMATRIX_H_
#define MYMATRIX_H_

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <limits>
#include <lapacke.h>
#include "myfun.h"

using namespace std;

#ifndef MyEPS_
#define MyEPS_
static double EPS = numeric_limits<double>::epsilon();
#endif
// =====================================================================
template<typename T>
struct brace {
	T first, second;
};

template<typename T> brace<T> make_brace(const T &m1, const T &m2) {
	brace<T> ans;
	ans.first = m1;
	ans.second = m2;
	return ans;
}

template<typename T>
struct triplet {
	T first, second, third;
};

template<typename T> triplet<T> make_triplet(const T &m1, const T &m2, 
		const T &m3) {
	triplet<T> ans;
	ans.first = m1;
	ans.second = m2;
	ans.third = m3;
	return ans;
}
// =====================================================================
class Matrix;
Matrix column_vector(const vector<double> &vec);
Matrix row_vector(const vector<double> &vec);
Matrix identityMatrix(const size_t &N);
Matrix diagonalMatrix(const vector<double> &vec);

class Matrix { // Matrix A_mXn
	friend Matrix operator+(const double &lhs, const Matrix &rhs); 
	friend Matrix operator-(const double &lhs, const Matrix &rhs);
	friend Matrix operator*(const double &lhs, const Matrix &rhs); 
	friend Matrix operator/(const double &lhs, const Matrix &rhs);
	friend ostream& operator<<(ostream &os, const Matrix& M);
  protected:
	vector<double> elements;
	void badConvergence(const string &str) const;
	void invalidArgument(const int &errid, const string &str) const;
	void checkDim(const size_t &i_, const string &str = "") const;
	void checkDim(const size_t &i_, const size_t &j_, 
		const string &str = "") const;
	void checkSize(const size_t &i_, const size_t &j_, 
		const string &str = "") const;
	void ensureSquare(const string &str = "") const;
	void ensureSameSize(const Matrix &B, const string &str = "") const;
	void ensureNonSingularity(const string &str = "") const;
	void ensureVector(const string &str) const;
  public:
	size_t m; // number of rows of the matrix
	size_t n; // number of columns of the matrix
	Matrix(); //constructor
    Matrix(const size_t &m_, const size_t &n_, const double &cns = 0.);
    Matrix(const size_t &m_, const size_t &n_, const double *elem, 
		const size_t &sz_elem);
	Matrix(const size_t &m_, const size_t &n_, 
		const vector<double> &elem);
    void swap(Matrix &);
    ~Matrix() {};//destructor
    void print() const;
    string str() const;
    double& at(const size_t &i_);
    const double& at(const size_t &i_) const;
    double& at(const size_t &i_, const size_t &j_);
    const double& at(const size_t &i_, const size_t &j_) const;
    vector<double> as_vector() const;
    vector<double> row(const size_t &r_) const;
    void assign2row(const size_t &r_, const double &value);
    void assign2row(const size_t &r_, const vector<double> &values);
    void assign2row(const size_t &r_, const double *values, 
		const size_t &sz);
    vector<double> column(const size_t &c_) const;
    void assign2column(const size_t &c_, const double &value);
    void assign2column(const size_t &c_, const vector<double> &values);
    void assign2column(const size_t &c_, const double *values, 
		const size_t &sz);
    vector<double> diagonal() const;
    void assign2diagonal(const double &value);
    void assign2diagonal(const vector<double> &values);
    void assign2diagonal(const double *values, const size_t &sz);
	void swap_rows(const size_t &r1, const size_t &r2);
	void swap_columns(const size_t &c1, const size_t &c2);

    void self_concat_row(const double &cte);
    void self_concat_row(const vector<double> &values);
	void self_concat_row(const Matrix &B);
	void self_concat_column(const double &cte);
	void self_concat_column(const vector<double> &values);
	void self_concat_column(const Matrix &B);
	void self_delete_row(const size_t &r_);
	void self_delete_column(const size_t &c_);
	
	Matrix concat_row(const double &cte) const;
    Matrix concat_row(const vector<double> &values) const;
	Matrix concat_row(const Matrix &B) const;
	Matrix concat_column(const double &cte) const;
	Matrix concat_column(const vector<double> &values) const;
	Matrix concat_column(const Matrix &B) const;
	Matrix delete_row(const size_t &r_) const;
	Matrix delete_column(const size_t &c_) const;

	size_t rank() const;
	double determinant() const;
	double trace() const;
	double norm() const;
	double infinity_norm() const;
	double frobenius_norm() const;
	double condition_number() const;
	double cofactor(const size_t &r_, const size_t &c_) const;
	double highest_eigenvalue() const;
    Matrix transpose() const;
    Matrix adjoint() const;
    Matrix inverse() const;
    brace<Matrix> eigen() const;
    triplet<Matrix> lu_p_factorization() const;
    brace<Matrix> qr_factorization() const;
    brace<Matrix> lq_factorization() const;
    triplet<Matrix> svd() const;
    Matrix pseudoinverse() const;
        
    Matrix operator+(const double &value) const;
    void operator+=(const double &value);
    Matrix operator-(const double &value) const;
    void operator-=(const double &value);
    Matrix operator*(const double &value) const;
    void operator*=(const double &value);
    Matrix operator/(const double &value) const;
    void operator/=(const double &value);
    Matrix operator+(const Matrix &value) const;
    void operator+=(const Matrix &value);
    Matrix operator-(const Matrix &value) const;
    void operator-=(const Matrix &value);
    Matrix operator*(const Matrix &value) const;
    void operator*=(const Matrix &value);
    bool operator==(const Matrix &B) const;
    bool operator!=(const Matrix &B) const;
    
    double min_all(const bool &absolute) const;
    double max_all(const bool &absolute) const;
    
    bool isFullRank() const;
    bool isOrthogonal() const;
    bool isDiagonal() const;
    bool isPositiveDefinite() const;
    bool isSingular() const;
    bool isSymmetric() const;
    bool isSkewSymmetric() const;
    bool isValid() const;
    
    void adjust_elements(const double &tol);
    
    void saveAs(const string &namefile) const;
};

Matrix::Matrix() : m(0), n(0) {
}

Matrix::Matrix(const size_t &m_, const size_t &n_, const double &cns) :
		m(m_), n(n_) {
	this->elements.assign(m*n,cns);
}

Matrix::Matrix(const size_t &m_, const size_t &n_, const double *elem, 
		const size_t &sz_elem) : m(m_), n(n_) {
	this->checkSize(this->m*this->n, sz_elem,"(CONSTRUCTOR)");
	this->elements.assign(elem, elem+sz_elem);
}

Matrix::Matrix(const size_t &m_, const size_t &n_, const vector<double> 
		&elem) : m(m_), n(n_) {
	this->checkSize(this->m*this->n, elem.size(),"(CONSTRUCTOR)");
	this->elements = elem;
}

void Matrix::swap(Matrix &other) {
    // need to swap all data members
    std::swap(this->elements, other.elements);
    std::swap(this->m, other.m);
    std::swap(this->n, other.n);
}
// =====================================================================
void Matrix::badConvergence(const string &str) const {
	string MSG = (str.empty() ? "" : "in ") + str + 
		(str.empty() ? "" : ", ") + 
		"the algorithm failed to converge.";
	throw std::runtime_error(MSG.c_str());
}

void Matrix::invalidArgument(const int &errid, const string &str) 
		const {
	string MSG = (str.empty() ? "" : "in ") + str + 
		(str.empty() ? "" : ", ") + 
		"invalid argument";
	cerr << "[ " << errid << " ] ";
	throw std::invalid_argument(MSG.c_str());
}

void Matrix::checkDim(const size_t &i_, const string &str) const {
	string MSG = (str.empty() ? "" : "in ") + str + 
		(str.empty() ? "" : ", ") + 
		"index exceeds matrix dimension";
	if(i_ >= this->elements.size()) 
		throw std::out_of_range(MSG.c_str());
}

void Matrix::checkDim(const size_t &i_, const size_t &j_, 
		const string &str) const {
	string MSG = (str.empty() ? "" : "in ") + str + 
		(str.empty() ? "" : ", ") + 
		"index exceeds matrix dimension";
	if(i_ >= this->m || j_ >= this->n) 
		throw std::out_of_range(MSG.c_str());
}

void Matrix::checkSize(const size_t &i_, const size_t &j_, 
		const string &str) 	const {
	string MSG = (str.empty() ? "" : "in ") + str + 
		(str.empty() ? "" : ", ") + 
		"sizes of matrices don't match";
	if(i_ != j_) 
		throw std::invalid_argument(MSG.c_str());
}

void Matrix::ensureSquare(const string &str) const {
	string MSG = (str.empty() ? "" : "in ") + str + 
		(str.empty() ? "" : ", ") + 
		"exclusive to square matrices";
	if(this->m != this->n) 
		throw std::invalid_argument(MSG.c_str());
}

void Matrix::ensureSameSize(const Matrix &B, const string &str) const {
	string MSG = (str.empty() ? "" : "in ") + str + 
		(str.empty() ? "" : ", ") + 
		"matrices must be of same size";
	if(this->m != B.m || this->n != B.n) 
		throw std::invalid_argument(MSG.c_str());
}

void Matrix::ensureNonSingularity(const string &str) const {
	string MSG = (str.empty() ? "" : "in ") + str + 
		(str.empty() ? "" : ", ") + 
		"matrix is singular";
	if(this->isSingular()) 
		throw std::invalid_argument(MSG.c_str());
}

void Matrix::ensureVector(const string &str) const {
	string MSG = (str.empty() ? "" : "in ") + str + 
		(str.empty() ? "" : ", ") + 
		"matrices must be row or column vectors";
	if(this->m != 1 && this->n != 1) 
		throw std::invalid_argument(MSG.c_str());
}
// =====================================================================
void Matrix::print() const {
	cout << endl << ">> Matrix_" << this->m << "X" << this->n 
		 << " =" << endl;
	for(size_t i = 0; i < this->m; i++) {
		for(size_t j = 0; j < this->n; j++)
			cout << this->elements.at(i*this->n + j) << '\t';
		cout << endl;
	}
	cout << endl;
}

string Matrix::str() const {
	string ST = "";
	for(size_t i = 0; i < this->m; i++) {
		for(size_t j = 0; j < this->n; j++)
			ST += my::num2str<double>(this->elements.at(i*this->n + j)) 
				+ ((j == this->n-1) ? '\0' : '\t');
		ST += '\n';
	}
	return ST;
}

ostream& operator<<(ostream &os, const Matrix& M) {
	return os << M.str();
}

double& Matrix::at(const size_t &i_) {
	this->checkDim(i_,"AT");
	return this->elements.at(i_);
}

const double& Matrix::at(const size_t &i_) const {
	this->checkDim(i_,"AT");
	return this->elements.at(i_);
}

double& Matrix::at(const size_t &i_, const size_t &j_) {
	this->checkDim(i_,j_,"AT");
	return this->elements.at(i_*(this->n)+j_);
}

const double& Matrix::at(const size_t &i_, const size_t &j_) const {
	this->checkDim(i_,j_,"AT");
	return this->elements.at(i_*(this->n)+j_);
}

vector<double> Matrix::as_vector() const {
	return this->elements;
}

vector<double> Matrix::row(const size_t &r_) const {
	this->checkDim(r_,0,"ROW");
	vector<double> row_(this->elements.begin()+(r_*this->n),
		this->elements.begin()+(r_*this->n+this->n));
	return row_;
}

void Matrix::assign2row(const size_t &r_, const double &value) {
	this->checkDim(r_,0,"ASSIGN2ROW");
	for(size_t i = 0; i < this->n; i++)
		this->elements.at(r_*this->n+i) = value;
}

void Matrix::assign2row(const size_t &r_, 
		const vector<double> &values) {
	this->checkDim(r_,0,"ASSIGN2ROW");
	this->checkSize(this->n,values.size(),"ASSIGN2ROW");
	for(size_t i = 0; i < this->n; i++)
		this->elements.at(r_*this->n+i) = values.at(i);
}

void Matrix::assign2row(const size_t &r_, const double *values, 
		const size_t &sz) {
	vector<double> vec(values,values+sz);
	this->assign2row(r_,vec);
}

vector<double> Matrix::column(const size_t &c_) const {
	this->checkDim(0,c_,"COLUMN");
	vector<double> column_(this->m,0.);
	for(size_t i = 0; i < this->m; i++)
		column_.at(i) = this->elements.at(i*this->n+c_);
	return column_;
}

void Matrix::assign2column(const size_t &c_, const double &value) {
	this->checkDim(0,c_,"ASSIGN2COLUMN");
	for(size_t i = 0; i < this->m; i++)
		this->elements.at(i*this->n+c_) = value;
}

void Matrix::assign2column(const size_t &c_, 
		const vector<double> &values) {
	this->checkDim(0,c_,"ASSIGN2COLUMN");
	this->checkSize(this->m,values.size(),"ASSIGN2COLUMN");
	for(size_t i = 0; i < this->m; i++)
		this->elements.at(i*this->n+c_) = values.at(i);
}

void Matrix::assign2column(const size_t &c_, const double *values, 
		const size_t &sz) {
	vector<double> vec(values,values+sz);
	this->assign2column(c_,vec);
}

vector<double> Matrix::diagonal() const {
	size_t dim = min(this->m,this->n);
	vector<double> diag_(dim,0.);
	for(size_t i = 0; i < dim; i++)
		diag_.at(i) = this->elements.at(i*this->n+i);
	return diag_;
}

void Matrix::assign2diagonal(const double &value) {
	size_t dim = min(this->m,this->n);
	for(size_t i = 0; i < dim; i++)
		this->elements.at(i*this->n+i) = value;
}

void Matrix::assign2diagonal(const vector<double> &values) {
	size_t dim = min(this->m,this->n);
	this->checkSize(dim,values.size(),"ASSIGN2DIAGONAL");
	for(size_t i = 0; i < dim; i++)
		this->elements.at(i*this->n+i) = values.at(i);
}

void Matrix::assign2diagonal(const double *values, const size_t &sz) {
	vector<double> vec(values,values+sz);
	this->assign2diagonal(vec);
}

void Matrix::swap_rows(const size_t &r1, const size_t &r2) {
	this->checkDim(r1,0,"SWAP_ROW");
	this->checkDim(r2,0,"SWAP_ROW");
	double aux;
	for(size_t j = 0; j < this->n; j++) {
		aux = this->elements.at(r1*this->n+j);
		this->elements.at(r1*this->n+j) = 
			this->elements.at(r2*this->n+j);
		this->elements.at(r2*this->n+j) = aux;
	}
}

void Matrix::swap_columns(const size_t &c1, const size_t &c2) {
	this->checkDim(0,c1,"SWAP_COLUMN");
	this->checkDim(0,c2,"SWAP_COLUMN");
	double aux;
	for(size_t i = 0; i < this->m; i++) {
		aux = this->elements.at(i*this->n+c1);
		this->elements.at(i*this->n+c1) = 
			this->elements.at(i*this->n+c2);
		this->elements.at(i*this->n+c2) = aux;
	}
}

void Matrix::self_concat_row(const double &cte) {
	Matrix M(this->concat_row(cte));
	this->swap(M);
}

void Matrix::self_concat_row(const vector<double> &values) {
	Matrix M(this->concat_row(values));
	this->swap(M);
}

void Matrix::self_concat_row(const Matrix &B) {
	Matrix M(this->concat_row(B));
	this->swap(M);
}

void Matrix::self_concat_column(const double &cte) {
	Matrix M(this->concat_column(cte));
	this->swap(M);
}

void Matrix::self_concat_column(const vector<double> &values) {
	Matrix M(this->concat_column(values));
	this->swap(M);
}

void Matrix::self_concat_column(const Matrix &B) {
	Matrix M(this->concat_column(B));
	this->swap(M);
}

void Matrix::self_delete_row(const size_t &r_) {
	Matrix M(this->delete_row(r_));
	this->swap(M);
}

void Matrix::self_delete_column(const size_t &c_) {
	Matrix M(this->delete_column(c_));
	this->swap(M);
}

Matrix Matrix::concat_row(const double &cte) const {
	vector<double> values(this->n,cte);
	return concat_row(values);
}

Matrix Matrix::concat_row(const vector<double> &values) const {
	this->checkSize(this->n,values.size(),"CONCAT_ROW");
	Matrix M(this->m+1,this->n);
	for(size_t i = 0; i < this->m; i++)
		for(size_t j = 0; j < M.n; j++)
			M.elements.at(i*M.n+j) = this->elements.at(i*this->n+j);
	for(size_t j = 0; j < M.n; j++)
		M.elements.at(this->m*M.n+j) = values.at(j);
	return M;
}

Matrix Matrix::concat_row(const Matrix &B) const {
	this->checkSize(this->n,B.n,"CONCAT_ROW");
	Matrix M(this->m+B.m,this->n);
	for(size_t i = 0; i < this->m; i++)
		for(size_t j = 0; j < M.n; j++)
			M.elements.at(i*M.n+j) = this->elements.at(i*this->n+j);
	for(size_t i = 0; i < B.m; i++)
		for(size_t j = 0; j < M.n; j++)
			M.elements.at((this->m+i)*M.n+j) = 
				B.elements.at(i*B.n+j);
	return M;
}

Matrix Matrix::concat_column(const double &cte) const {
	vector<double> values(this->m,cte);
	return concat_column(values);
}

Matrix Matrix::concat_column(const vector<double> &values) const {
	this->checkSize(this->m,values.size(),"CONCAT_COLUMN");
	Matrix M(this->m,this->n+1);
	for(size_t i = 0; i < M.m; i++)
		for(size_t j = 0; j < this->n; j++)
			M.elements.at(i*M.n+j) = this->elements.at(i*this->n+j);
	for(size_t i = 0; i < M.m; i++)
		M.elements.at(i*M.n+this->n) = values.at(i);
	return M;
}

Matrix Matrix::concat_column(const Matrix &B) const {
	this->checkSize(this->m,B.m,"CONCAT_COLUMN");
	Matrix M(this->m,this->n+B.n);
	for(size_t i = 0; i < M.m; i++)
		for(size_t j = 0; j < this->n; j++)
			M.elements.at(i*M.n+j) = this->elements.at(i*this->n+j);
	for(size_t i = 0; i < M.m; i++)
		for(size_t j = 0; j < B.n; j++)
			M.elements.at(i*M.n+(this->n+j)) = B.elements.at(i*B.n+j);
	return M;
}

Matrix Matrix::delete_row(const size_t &r_) const {
	this->checkDim(r_,0,"DELETE_ROW");
	Matrix M(this->m-1,this->n);
	for(size_t i = 0; i < this->m; i++)
		for(size_t j = 0; j < this->n; j++)
			if(i != r_)
				M.elements.at(((i<r_)?i:i-1)*M.n+j) = 
					this->elements.at(i*this->n+j);
	return M;
}

Matrix Matrix::delete_column(const size_t &c_) const {
	this->checkDim(0,c_,"DELETE_COLUMN");
	Matrix M(this->m,this->n-1);
	for(size_t i = 0; i < this->m; i++)
		for(size_t j = 0; j < this->n; j++)
			if(j != c_)
				M.elements.at(i*M.n+((j<c_)?j:j-1)) = 
					this->elements.at(i*this->n+j);
	return M;
}
// ===================================
size_t Matrix::rank() const {
	size_t rk = 0;
	this->checkSize((this->m)*(this->n),this->elements.size(),"RANK");
	lapack_int m = this->m, n = this->n;
	lapack_int lda = n, ldu = m, ldvt = n, info;
	double *superb = new double[min(m,n)-1];
	double *a = new double[m*n],
		   *s = new double[min(m,n)], 
		   *u = new double[m*m], 
		   *vt = new double[n*n];
	for(size_t i = 0; i < this->elements.size(); i ++)
		a[i] = this->elements.at(i);
	info = LAPACKE_dgesvd( LAPACK_ROW_MAJOR, 'N', 'N', m, n, a, lda,
		s, u, ldu, vt, ldvt, superb );
    if(info < 0) 	// Check for invalid (NaN)
		this->invalidArgument((int)info,"RANK");
	if(info > 0) 	// Check for convergence
		this->badConvergence("RANK");
	for(int i = 0; i < min(m,n); i++)
		rk += (my::isLikelyZero(s[i],max(this->m,this->n))) ? 0 : 1;
	delete [] superb;
	delete [] a;
	delete [] s;
	delete [] u;
	delete [] vt;
	return rk;
}

double Matrix::determinant() const {
	this->ensureSquare("DETERMINANT");
	double det;
	switch(this->m) {
		case 1:
			det = this->elements.at(0);
			break;
		case 2:
			det = this->elements.at(0)*this->elements.at(3) -
				this->elements.at(1)*this->elements.at(2);
			break;
		case 3:
			det = this->elements.at(0)*this->elements.at(4)*
					this->elements.at(8) 
				+ this->elements.at(1)*this->elements.at(5)*
					this->elements.at(6) 
				+ this->elements.at(2)*this->elements.at(3)*
					this->elements.at(7) 
				- this->elements.at(2)*this->elements.at(4)*
					this->elements.at(6)
				- this->elements.at(0)*this->elements.at(5)*
					this->elements.at(7)
				- this->elements.at(1)*this->elements.at(3)*
					this->elements.at(8);
			break;
		default:
			det = 0;
			for(size_t j = 0; j < this->n; j++)
				det += this->cofactor(0,j)*this->elements.at(j); // i=0 
	}
	return det;
}

double Matrix::trace() const {
	this->ensureSquare("TRACE");
	vector<double> diag(this->diagonal());
	double total(0.);
	for(vector<double>::iterator it = diag.begin(); 
		it < diag.end(); it++) total += *it;
	return total;
}

double Matrix::norm() const {
	// for square matrices, this is the same as the "spectral norm",
	// the "induced matrix 2-norm" or the "Euclidean norm"
	if(this->m != 1 && this->n != 1) {
		this->checkSize((this->m)*(this->n),this->elements.size(),
			"NORM");
		lapack_int m = this->m, n = this->n;
		lapack_int lda = n, ldu = m, ldvt = n, info;
		double *superb = new double[min(m,n)-1];
		double *a = new double[m*n],
			   *s = new double[min(m,n)], 
			   *u = new double[m*m], 
			   *vt = new double[n*n];
		for(size_t i = 0; i < this->elements.size(); i ++)
			a[i] = this->elements.at(i);
		info = LAPACKE_dgesvd( LAPACK_ROW_MAJOR, 'N', 'N', m, n, a, lda,
			s, u, ldu, vt, ldvt, superb );
		if(info < 0) 	// Check for invalid (NaN)
			this->invalidArgument((int)info,"NORM");
		if(info > 0) 	// Check for convergence
			this->badConvergence("NORM");
		double hSingVal = s[0];
		for(int i = 1; i < min(m,n); i++)
			if(s[i] > hSingVal)
				hSingVal = s[i];
		delete [] superb;
		delete [] a;
		delete [] s;
		delete [] u;
		delete [] vt;
		return hSingVal;
	} else
		return sqrt(((this->m == 1) ? (*this)*(this->transpose()) :
			(this->transpose())*(*this)).elements.at(0));
}

double Matrix::infinity_norm() const {
	// Maximum absolute row sum norm
	double infnorm = 0., sum_;
	for(size_t i = 0; i < this->m; i++) {
		sum_ = 0.;
		for(size_t j = 0; j < this->n; j++)
			sum_ += abs(this->elements.at(i*this->n+j));
		if(sum_ > infnorm)
			infnorm = sum_;
	}
	return infnorm;
}

double Matrix::frobenius_norm() const {
	this->checkSize((this->m)*(this->n),this->elements.size(),
		"FRBNORM");
	double sum = 0.;
	for(size_t i = 0; i < this->m; i++)
		for(size_t j = 0; j < this->n; j++)
			sum += pow(this->elements.at(i*this->n+j),2.);
	return sqrt(sum);
}

double Matrix::condition_number() const
{
	triplet<Matrix> SVD = this->svd();
	vector<double> singular_values = SVD.second.diagonal();
	return my::max(singular_values)/my::min(singular_values);
}
double Matrix::cofactor(const size_t &r_, const size_t &c_) const {
	this->ensureSquare("COFACTOR");
	Matrix CF(*this);
	CF.self_delete_row(r_);
	CF.self_delete_column(c_);
	return (((r_+c_)%2)?-1:1)*CF.determinant();
}

Matrix Matrix::transpose() const {
	Matrix T(this->n,this->m);
	for(size_t i = 0; i < this->m; i++)
		for(size_t j = 0; j < this->n; j++)
			T.elements.at(j*this->m+i) = this->elements.at(i*this->n+j);
	return T;
}

Matrix Matrix::adjoint() const {
	this->ensureSquare("ADJOINT");
	Matrix A(this->m,this->n);
	for(size_t i = 0; i < this->m; i++)
		for(size_t j = 0; j < this->n; j++)
			A.elements.at(i*this->n+j) = cofactor(i,j);
	return A.transpose();
}

Matrix Matrix::inverse() const {
	this->ensureSquare("INVERSE");
	this->ensureNonSingularity("INVERSE");
	return this->adjoint()/this->determinant();
}

brace<Matrix> Matrix::eigen() const {
	this->checkSize((this->n)*(this->n),this->elements.size(),"EIGEN");
	this->ensureSquare("EIGEN");
	lapack_int n = this->n, info;
	double *a = new double[n*n],
		   *b = new double[n*n],
		   *alphar = new double[n],
		   *alphai = new double[n],
		   *beta = new double[n],
		   *vl = new double,
		   *vr = new double[n*n];
	for(size_t i = 0; i < this->elements.size(); i ++) {
		a[i] = this->elements.at(i);
		b[i] = (i%(n+1) == 0) ? 1. : 0.;
	}
	info = LAPACKE_dggev( LAPACK_ROW_MAJOR, 'N', 'V', n, a, n, b, n, 
		alphar, alphai, beta, vl, 1, vr, n );
	if(info < 0) 	// Check for invalid (NaN)
		this->invalidArgument((int)info,"EIGEN");
	if(info > 0) 	// Check for convergence
		this->badConvergence("EIGEN");
	Matrix EVAL(n,n),
		   EVEC(n,n,vr,n*n);
	vector<double> lambda(n,0.);
	for(int i = 0; i < n; i++) {
		lambda.at(i) = (abs(alphai[i]) <= EPS) ? ((abs(beta[i]) <= EPS) 
			? NAN : alphar[i]/beta[i] ) : -NAN;
	}
	EVAL.assign2diagonal(lambda);
	delete [] a;
	delete [] b;
	delete [] alphar;
	delete [] alphai;
	delete [] beta;
	delete vl;
	delete [] vr;
	return make_brace<Matrix>(EVAL,EVEC);
}

triplet<Matrix> Matrix::lu_p_factorization() const {
	lapack_int m = this->m, n = this->n, info;
	double *a = new double[m*n];
	lapack_int *ipiv = new lapack_int[max(1,min(m,n))];
	for(size_t i = 0; i < this->elements.size(); i ++) {
		a[i] = this->elements.at(i);
	}
	info = LAPACKE_dgetrf( LAPACK_ROW_MAJOR, m, n, a, n, ipiv );
	if(info < 0) 	// Check for invalid (NaN)
		this->invalidArgument((int)info,"LU_P_FACTORIZATION");
	if(info > 0) 	// Check for convergence
		cerr << ">> Warning! Matrix U is singular." << endl;
	Matrix L(m,min(m,n)), 
		   U(min(m,n),n), 
		   P(m,m);
	vector<double> ones(min(m,n),1.);
	L.assign2diagonal(ones);
	for(int i = 0; i < m; i++)
		for(int j = 0; j < n; j++)
			if(j >= i)
				U.at(i,j) = a[i*n+j];
			else
				L.at(i,j) = a[i*n+j];
	vector<size_t> seq;
	size_t aux;
	for(int i = 0; i < m; i++)
		seq.push_back(i);
	for(int i = 0; i < min(m,n); i++) {
		aux = seq.at(i);
		seq.at(i) = seq.at(ipiv[i]-1);
		seq.at(ipiv[i]-1) = aux;
	}
	for(int i = 0; i < m; i++)
		P.at(i,seq.at(i)) = 1;
	delete [] a;
	delete [] ipiv;
	return make_triplet<Matrix>(L,U,P);
}

brace<Matrix> Matrix::qr_factorization() const {
	lapack_int m = this->m, n = this->n, info;
	double *a = new double[m*n];
	double *tau = new double[max(1,min(m,n))];
	for(size_t i = 0; i < this->elements.size(); i ++) {
		a[i] = this->elements.at(i);
	}
    info = LAPACKE_dgeqrf( LAPACK_ROW_MAJOR, m, n, a, n, tau );
	if(info < 0) 	// Check for invalid (NaN)
		this->invalidArgument((int)info,"QR_FACTORIZATION");
	/* The matrix Q is represented as a product of elementary reflectors
	*
	*     Q = H(1) H(2) . . . H(k), where k = min(m,n).
	*
	*  Each H(i) has the form
	*
	*     H(i) = I - tau * v * v'
	*
	*  where tau is a real scalar, and v is a real vector with
	*  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in 
	*  A(i+1:m,i), and tau in TAU(i).*/	
	Matrix Q(m,min(m,n)), 
		   R(min(m,n),n), 
		   H(identityMatrix(m)),
		   v(m,1);
	vector<double> ones(min(m,n),1.);
	Q.assign2diagonal(ones);
	for(int i = 0; i < m; i++)
		for(int j = 0; j < n; j++)
			if(j >= i)
				R.at(i,j) = a[i*n+j];
			else
				Q.at(i,j) = a[i*n+j];
	for(int i = 0; i < min(m,n); i++) {
		v = column_vector(Q.column(i));
		H *= identityMatrix(m)-(v*v.transpose())*tau[i];
	}
	Q = H;
	if(m > n)
		for(int i = n; i < m; i++)
			Q.self_delete_column(Q.n-1);
	delete [] a;
	delete [] tau;
	return make_brace<Matrix>(Q,R);
}

brace<Matrix> Matrix::lq_factorization() const {
	lapack_int m = this->m, n = this->n, info;
	double *a = new double[m*n];
	double *tau = new double[max(1,min(m,n))];
	for(size_t i = 0; i < this->elements.size(); i ++) {
		a[i] = this->elements.at(i);
	}
    info = LAPACKE_dgelqf( LAPACK_ROW_MAJOR, m, n, a, n, tau );
	if(info < 0) 	// Check for invalid (NaN)
		this->invalidArgument((int)info,"LQ_FACTORIZATION");
	Matrix L(m,min(m,n)), 
		   Q(min(m,n),n), 
		   H(identityMatrix(n)),
		   v(n,1);
	vector<double> ones(min(m,n),1.);
	Q.assign2diagonal(ones);
	for(int i = 0; i < m; i++)
		for(int j = 0; j < n; j++)
			if(j <= i)
				L.at(i,j) = a[i*n+j];
			else
				Q.at(i,j) = a[i*n+j];
	for(int i = 0; i < min(m,n); i++) {
		v = row_vector(Q.row(i));
		H *= identityMatrix(n)-(v.transpose()*v)*tau[i];
	}
	Q = H.transpose();
	if(m < n)
		for(int i = m; i < n; i++)
			Q.self_delete_row(Q.m-1);
	delete [] a;
	delete [] tau;
	return make_brace<Matrix>(L,Q);
}

triplet<Matrix> Matrix::svd() const {
	this->checkSize((this->m)*(this->n),this->elements.size(),"SVD");
	lapack_int m = this->m, n = this->n;
	lapack_int lda = n, ldu = m, ldvt = n, info;
	double *superb = new double[min(m,n)-1];
	double *a = new double[m*n],
		   *s = new double[min(m,n)], 
		   *u = new double[m*m], 
		   *vt = new double[n*n];
	for(size_t i = 0; i < this->elements.size(); i ++)
		a[i] = this->elements.at(i);
	info = LAPACKE_dgesvd( LAPACK_ROW_MAJOR, 'A', 'A', m, n, a, lda,
		s, u, ldu, vt, ldvt, superb );
    if(info < 0) 	// Check for invalid (NaN)
		this->invalidArgument((int)info,"SVD");
	if(info > 0) 	// Check for convergence
		this->badConvergence("SVD");
	// Adjust zeros
	for(size_t i = 0; i < min(this->m,this->n); i++)
		if(my::isLikelyZero(s[i],max(this->m,this->n)))
			s[i] = 0.;
	for(int i = 0; i < m*m; i++)
		if(my::isLikelyZero(u[i],max(this->m,this->n)))
			u[i] = 0.;
	for(int i = 0; i < n*n; i++)
		if(my::isLikelyZero(vt[i],max(this->m,this->n)))
			vt[i] = 0.;
	// Create matrices
	Matrix S(m,n),
		   U(m,m,u,m*m),
		   VT(n,n,vt,n*n);
	S.assign2diagonal(s,min(m,n));
	delete [] superb;
	delete [] a;
	delete [] s;
	delete [] u;
	delete [] vt;
	// Return matrices U, Sigma and V
	// To retrieve original matrix: U*Sigma*V.transpose()
	return make_triplet<Matrix>(U,S,VT.transpose());
}

Matrix Matrix::pseudoinverse() const {
	triplet<Matrix> SVD = this->svd();
	vector<double> sv = SVD.second.diagonal();
	double maxsv = 0;
	for(size_t i = 0; i < sv.size(); i++)
		if(maxsv < sv.at(i))
			maxsv = sv.at(i);
	for(size_t i = 0; i < sv.size(); i++) {
		if(sv.at(i) <= this->m*this->n*EPS)
			sv.at(i) = 0.;
		else
			sv.at(i) = 1./sv.at(i);
	}
	SVD.second.assign2diagonal(sv);
	return (SVD.first*SVD.second*SVD.third.transpose()).transpose();
}
// =========
double Matrix::min_all(const bool &absolute) const
{
	double minelem = (absolute)? abs(this->elements.front()) : 
		this->elements.front();
	for(size_t i = 1; i < this->elements.size(); i++)
		if((absolute)? abs(this->elements.at(i)) < minelem : 
			this->elements.at(i) < minelem)
				minelem = (absolute)? abs(this->elements.at(i)) : 
					this->elements.at(i);
	return minelem;
}
double Matrix::max_all(const bool &absolute) const
{
	double maxelem = (absolute)? abs(this->elements.front()) : 
		this->elements.front();
	for(size_t i = 1; i < this->elements.size(); i++)
		if((absolute)? abs(this->elements.at(i)) > maxelem : 
			this->elements.at(i) > maxelem)
				maxelem = (absolute)? abs(this->elements.at(i)) : 
					this->elements.at(i);
	return maxelem;
}
// =========
Matrix Matrix::operator+(const double &value) const {
	Matrix A(*this);
	A += value;
	return A;
}

void Matrix::operator+=(const double &value) {
	for(size_t i = 0; i < this->m; i++)
		for(size_t j = 0; j < this->n; j++)
			this->elements.at(i*this->n+j) += value;
}

Matrix Matrix::operator-(const double &value) const {
	Matrix A(*this);
	A -= value;
	return A;
}

void Matrix::operator-=(const double &value) {
	for(size_t i = 0; i < this->m; i++)
		for(size_t j = 0; j < this->n; j++)
			this->elements.at(i*this->n+j) -= value;
}

Matrix Matrix::operator*(const double &value) const {
	Matrix A(*this);
	A *= value;
	return A;
}

void Matrix::operator*=(const double &value) {
	for(size_t i = 0; i < this->m; i++)
		for(size_t j = 0; j < this->n; j++)
			this->elements.at(i*this->n+j) *= value;
}

Matrix Matrix::operator/(const double &value) const {
	Matrix A(*this);
	A /= value;
	return A;
}

void Matrix::operator/=(const double &value) {
	for(size_t i = 0; i < this->m; i++)
		for(size_t j = 0; j < this->n; j++)
			this->elements.at(i*this->n+j) /= value;
}

Matrix Matrix::operator+(const Matrix &B) const {
	Matrix A(*this);
	A += B;
	return A;
}

void Matrix::operator+=(const Matrix &B) {
	this->ensureSameSize(B,"OPERATOR+");
	for(size_t i = 0; i < this->m; i++)
		for(size_t j = 0; j < this->n; j++)
			this->elements.at(i*this->n+j) += B.elements.
				at(i*this->n+j);
}

Matrix Matrix::operator-(const Matrix &B) const {
	Matrix A(*this);
	A -= B;
	return A;
}

void Matrix::operator-=(const Matrix &B) {
	this->ensureSameSize(B,"OPERATOR-");
	for(size_t i = 0; i < this->m; i++)
		for(size_t j = 0; j < this->n; j++)
			this->elements.at(i*this->n+j) -= B.elements.
				at(i*this->n+j);
}

Matrix Matrix::operator*(const Matrix &B) const {
	this->checkSize(this->n,B.m,"OPERATOR*");
	Matrix R(this->m,B.n);
	for(size_t i = 0; i < this->m; i++)
		for(size_t j = 0; j < this->n; j++)
			for(size_t k = 0; k < B.n; k++) {
				R.elements.at(i*B.n+k) += this->elements.
					at(i*this->n+j)*B.elements.at(j*B.n+k);
				if(abs(R.elements.at(i*B.n+k)) <= this->m*EPS) 
					R.elements.at(i*B.n+k) = 0.;
			}
	return R;
}

void Matrix::operator*=(const Matrix &B) {
	Matrix A = (*this)*B;
	this->swap(A);
}

bool Matrix::operator==(const Matrix &B) const {
	if(this->m != B.m || this->n != B.n)
		return false;
	else {
		bool chk = true;
		for(size_t i = 0; i < this->m; i++) {
			for(size_t j = 0; j < this->n; j++) {
				if(abs(this->elements.at(i*this->n+j)-
						B.elements.at(i*B.n+j)) > 0.5*EPS) {
					chk = false;
					break;
				}
			}
			if(!chk) break;
		}
		return chk;
	}
}

bool Matrix::operator!=(const Matrix &B) const {
	return !((*this) == B);
}
// ======
Matrix operator+(const double &lhs, const Matrix &rhs) {
	return rhs+lhs;
}
Matrix operator-(const double &lhs, const Matrix &rhs) {
	return rhs*(-1.)+lhs;
}
Matrix operator*(const double &lhs, const Matrix &rhs) {
	return rhs*lhs;
}
Matrix operator/(const double &lhs, const Matrix &rhs) {
	// This is not an inverse! It's an element-wise division
	Matrix result(rhs);
	for(size_t i = 0; i < result.m; i++)
		for(size_t j = 0; j < result.n; j++)
			result.at(i,j) = lhs/result.at(i,j);
	return result;
}
// ======
bool Matrix::isFullRank() const {
	return this->rank() == min(this->m,this->n);
}

bool Matrix::isOrthogonal() const {
	Matrix O = (*this)*(this->transpose());
	return (O == identityMatrix(O.m));
}

bool Matrix::isDiagonal() const {
	Matrix D = diagonalMatrix(this->diagonal());
	if(this->m != this->n) {
		if(this->m > this->n) {
			vector<double> zeros(this->n,0.);
			for(size_t i = this->n; i < this->m; i++)
				D.self_concat_row(zeros);
		} else {
			vector<double> zeros(this->m,0.);
			for(size_t i = this->m; i < this->n; i++)
				D.self_concat_column(zeros);
		}
	}
	return (*this == D);
}

bool Matrix::isPositiveDefinite() const {
	//real symmetric matrix
	if(!(this->isSymmetric()))
		return false;
	else {
		// if all its eigenvalues are positive.
		cout << "not implemented yet!" << endl;
		return false;
	}
}

bool Matrix::isSingular() const {
	return (my::isLikelyZero(this->determinant()));
}

bool Matrix::isSymmetric() const {
	return (*this == this->transpose());
}

bool Matrix::isSkewSymmetric() const {
	return (*this == this->transpose()*(-1.));
}

bool Matrix::isValid() const {
	bool chk = true;
	for(vector<double>::const_iterator it = this->elements.begin();
		it != this->elements.end(); it++)
	{
		if(std::isnan(*it) || my::isinf(*it))
		{
			chk = false;
			break;
		}
	}
	return chk;
}

void Matrix::adjust_elements(const double &tol)
{
	for(size_t i = 0; i < this->elements.size(); i++)
		if(my::isLikelyZero(this->elements.at(i),this->m*this->n) || 
				abs(this->elements.at(i)) < tol) 
			this->elements.at(i) = 0.;
}

void Matrix::saveAs(const string &namefile) const
{
	ofstream matrixfile;
	matrixfile.open(namefile);
	if (matrixfile.is_open())
	{
		for(size_t i = 0; i < this->m; i++)
		{
			for(size_t j = 0; j < this->n; j++)
				matrixfile << this->elements.at(i*(this->n)+j) 
						   << ((j+1 == this->n) ? "" : "\t");
			matrixfile << endl;
		}
	}
	else
	{
		cerr << "Unable to open " << namefile << " file!" << endl;
	}
	matrixfile.close();
}
// =====================================================================
Matrix column_vector(const vector<double> &vec) {
	Matrix A(vec.size(),1,vec);
	return A;
}

Matrix column_vector(const size_t &qty, const double &val) {
	vector<double> vec(qty,val);
	Matrix A(qty,1,vec);
	return A;
}

Matrix row_vector(const vector<double> &vec) {
	Matrix A(1,vec.size(),vec);
	return A;
}

Matrix row_vector(const size_t &qty, const double &val) {
	vector<double> vec(qty,val);
	Matrix A(1,qty,vec);
	return A;
}

Matrix subMatrix(const Matrix &M, const size_t &fstr, 
		const size_t &fstc, const size_t &lstr, const size_t &lstc) {
	cout << "row: " << fstr << " " << lstr << endl;
	cout << "col: " << fstc << " " << lstc << endl;
	cout << "M: " << M.m << " " << M.n << endl;
	if(fstr > lstr || lstr >= M.m || fstc > lstc || lstc >= M.n) {
		string MSG = "in SUBMATRIX, index exceeds matrix dimension";
		throw std::out_of_range(MSG.c_str());
	}
	Matrix A(lstr-fstr+1,lstc-fstc+1);
	cout << A.m << " " << A.n << endl;
	for(size_t i = 0; i < A.m; i++)
		for(size_t j = 0; j < A.n; j++)
			A.at(i,j) = M.at(i+fstr,j+fstc);
	return A;
}

Matrix identityMatrix(const size_t &N) {
	Matrix A(N,N);
	A.assign2diagonal(1.);
	return A;
}

Matrix diagonalMatrix(const vector<double> &vec) {
	Matrix A(vec.size(),vec.size());
	A.assign2diagonal(vec);
	return A;
}

Matrix linear_AASTtoRowVector_VAR(const AAST &tree, const size_t &n_var) 
{
	Matrix A, b(n_var, 1);
	vector<double> aux(n_var, 0.);
	double cnst = tree.eval_with_variables(aux);
	aux.at(0) = 1.;
	A = row_vector(aux);
	b.at(0) = tree.eval_with_variables(aux) - cnst;
	for(size_t i = 1; i < n_var; i++)
	{
		aux.assign(n_var, 0.);
		aux.at(i) = 1.;
		A.self_concat_row(row_vector(aux));
		b.at(i) = tree.eval_with_variables(aux) - cnst;
	}
	return ((A.pseudoinverse()*b).transpose()).concat_column(cnst);
}

Matrix linear_AASTtoRowVector_UNK(const AAST &tree, const size_t &n_unk) 
{
	Matrix A, b(n_unk, 1);
	vector<double> aux(n_unk, 0.);
	double cnst = tree.eval_with_unknowns(aux);
	aux.at(0) = 1.;
	A = row_vector(aux);
	b.at(0) = tree.eval_with_unknowns(aux) - cnst;
	for(size_t i = 1; i < n_unk; i++)
	{
		aux.assign(n_unk, 0.);
		aux.at(i) = 1.;
		A.self_concat_row(row_vector(aux));
		b.at(i) = tree.eval_with_unknowns(aux) - cnst;
	}
	return ((A.pseudoinverse()*b).transpose()).concat_column(cnst);
}
// =====================================================================
#endif
