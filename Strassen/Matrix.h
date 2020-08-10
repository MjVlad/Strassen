#ifndef _my_lab_h_
#define _my_lab_h_
#include<vector>
using namespace std;
class Matrix {
private:
	vector<vector<double>> mat1;
public:
	Matrix();
	Matrix(unsigned int a);
	void set(double val = 1);
	void print();
	int stroke();
	void knife();
	double det();
	void operator+=(Matrix right);
	void operator-=(Matrix right);
	int add_zero();
	void cut(Matrix& m11, Matrix& m12, Matrix& m21, Matrix& m22);
	void link(Matrix& m11, Matrix& m12, Matrix& m21, Matrix& m22);
	Matrix fin_cut(int size);
	friend Matrix operator*(Matrix& l, Matrix& r);
	friend Matrix multiStrassen(Matrix l, Matrix r, int mlt_thread);
};

Matrix operator+(Matrix l, Matrix r);
Matrix operator-(Matrix l, Matrix& r);


#endif // _my_lab_h_


