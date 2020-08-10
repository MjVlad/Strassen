#include<iostream>
#include<vector>
#include<cmath>
#include"Matrix.h"

using namespace std;

Matrix::Matrix(unsigned int a) {
	if (a < 1) {
		cout << endl << "You fucking idiot!!!";
		exit(-1);
	}
	mat1.resize(a);
	for (unsigned int i = 0; i < a; i++) {
		mat1.at(i).resize(a);
	}
}

void Matrix::set() {
	double tmp;
	double j = 1;
	for (unsigned int i = 0; i < mat1.size(); i++) {
		for (unsigned int k = 0; k < mat1.at(0).size(); k++) {
			//cout << "input [" << i + 1 << "] [" << k + 1 << "] ==>  ";
			//cin >> tmp;
			tmp = j*j;
			j++;
			mat1.at(i).at(k) = 1;
		}
	}
	cout << "\n=======================================================\n" << endl;
}

void Matrix::print() {
	for (unsigned int i = 0; i < mat1.size(); i++) {
		for (unsigned int k = 0; k < mat1.at(0).size(); k++) {
			cout /*<< "[" << i + 1 << "][" << k + 1 << "]<==>"*/ << mat1.at(i).at(k) << "    ";
		}
		cout << endl;
	}
	cout << "\n=======================================================\n" << endl;
}

int Matrix::stroke() {
	int size = mat1.at(0).size();
	if (mat1.at(0).at(0) == 0) {
		int i = 1;
		for (i; i < size; i++) {
			if (mat1.at(i).at(0) == 0)
				continue;
			else 
				break;
		}
		if (i == size) return 0;
		for (int j = 0; j < size; j++) {
			mat1.at(0).at(j) += mat1.at(i).at(j);
		}
	}
	for (int i = 1; i < size; i++) {
		double coe = mat1.at(i).at(0) / mat1.at(0).at(0);
		mat1.at(i).at(0) = 0;
		for (int j = 1; j < size; j++) {
			mat1.at(i).at(j) -= mat1[0][j] * coe;
		}
	}
	return 1;
}

void Matrix::knife() {
	int size = mat1.at(0).size();
	vector<vector<double>> mat2;
	mat2.resize(size - 1);
	for (int i = 1; i < size; i++) {
		mat2.at(i - 1).resize((size - 1));
		for (int j = 1; j < size; j++) {
			mat2.at(i - 1).at(j - 1) = mat1.at(i).at(j);
		}
		mat1.at(i).clear();
	}
	mat1 = mat2;
}

double Matrix::det() {
	int size = mat1.at(0).size();
	double det = 1;
	for (int i = 0; i < size - 2; i++) {
		cout << endl << i;
		int tmp = stroke();
		if (tmp == 0) return 0;
		det *= mat1.at(0).at(0);
		cout << "<<---->>" << det;
		knife();
	}
	det *= mat1.at(0).at(0) * mat1.at(1).at(1) - mat1.at(1).at(0) * mat1.at(0).at(1);
	return det;
}

void Matrix::operator+=(Matrix right) {
	int size = mat1.at(0).size();
	if (size != right.mat1.at(0).size()) return;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			mat1.at(i).at(j) += right.mat1.at(i).at(j);
		}
	}
	return;
}

void Matrix::operator-=(Matrix right) {
	int size = mat1.at(0).size();
	if (size != right.mat1.at(0).size()) return;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			mat1.at(i).at(j) -= right.mat1.at(i).at(j);
		}
	}
	return;
}

int Matrix::add_zero() {
	int size = 1;
	while (size < mat1.at(0).size()) {
		size <<= 1;
	}
	mat1.resize(size);
	for (int i = 0; i < size; i++) {
		mat1.at(i).resize(size);
	}
	return size;
}

void Matrix::cut(Matrix& m11, Matrix& m12, Matrix& m21, Matrix& m22) {
	int size = m11.mat1.at(0).size();
	int tmp = size * 2;
	for (int i = 0; i < size; i++) 
		for (int j = 0; j < size; j++) 
			m11.mat1.at(i).at(j) = mat1.at(i).at(j);
	for (int i = 0; i < size; i++) 
		for (int j = size; j < tmp; j++) 
			m12.mat1.at(i).at(j - size) = mat1.at(i).at(j);
	for (int i = size; i < tmp; i++) 
		for (int j = 0; j < size; j++) 
			m21.mat1.at(i - size).at(j) = mat1.at(i).at(j);
	for (int i = size; i < tmp; i++) 
		for (int j = size; j < tmp; j++) 
			m22.mat1.at(i - size).at(j - size) = mat1.at(i).at(j);
}

void Matrix::link(Matrix& m11, Matrix& m12, Matrix& m21, Matrix& m22) {
	int size = m11.mat1.at(0).size();
	int tmp = size * 2;
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			mat1.at(i).at(j) = m11.mat1.at(i).at(j);
	for (int i = 0; i < size; i++)
		for (int j = size; j < tmp; j++)
			mat1.at(i).at(j) = m12.mat1.at(i).at(j - size);
	for (int i = size; i < tmp; i++)
		for (int j = 0; j < size; j++)
			mat1.at(i).at(j) = m21.mat1.at(i - size).at(j);
	for (int i = size; i < tmp; i++)
		for (int j = size; j < tmp; j++)
			mat1.at(i).at(j) = m22.mat1.at(i - size).at(j - size);
}

Matrix Matrix::fin_cut(int size) {
	Matrix copy(size);
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			copy.mat1.at(i).at(j) = mat1.at(i).at(j);
		}
	}
	return copy;
}

Matrix operator*(Matrix& l, Matrix& r) {
	int size = l.mat1.at(0).size();
	Matrix m(size);
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			for (int k = 0; k < size; k++)
				m.mat1.at(i).at(j) += l.mat1.at(i).at(k) * r.mat1.at(k).at(j);
	return m;
}

Matrix multiStrassen(Matrix l, Matrix r) {
	static int size = 0;
	int tmp = l.mat1.at(0).size();
	if (size == 0) {
		size = l.add_zero();
		r.add_zero();
	}
	if (size <= 64) return l * r;
	Matrix out(size);
	size >>= 1;
	Matrix a11(size);
	Matrix a12(size);
	Matrix a21(size);
	Matrix a22(size);
	Matrix b11(size);
	Matrix b12(size);
	Matrix b21(size);
	Matrix b22(size);
	l.cut(a11, a12, a21, a22);
	r.cut(b11, b12, b21, b22);
	Matrix p1 = multiStrassen(a11, b12 - b22);
	Matrix p2 = multiStrassen(a11 + a12, b22);
	Matrix p3 = multiStrassen(a21 + a22, b11);
	Matrix p4 = multiStrassen(a22, b21 - b11);
	Matrix p5 = multiStrassen(a11 + a22, b11 + b22);
	Matrix p6 = multiStrassen(a12 - a22, b21 + b22);
	Matrix p7 = multiStrassen(a11 - a21, b11 + b12);
	Matrix c11 = (p5 + p4) + (p6 - p2);
	Matrix c12 = (p1 + p2);
	Matrix c21 = (p3 + p4);
	Matrix c22 = (p1 - p3) + (p5 - p7);
	out.link(c11, c12, c21, c22);
	out = out.fin_cut(tmp);
	return out;
}

Matrix operator+(Matrix l, Matrix r) {
	l += r;
	return l;
}

Matrix operator-(Matrix l, Matrix& r) {
	l -= r;
	return l;
}