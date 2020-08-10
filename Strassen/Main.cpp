#include "Matrix.h"
#include<iostream>
#include<vector>
#include <ctime>

using namespace std;

int main() {
	
	int size = 0;
	cin >> size;
	Matrix mat1 = Matrix(size);
	mat1.set(); 
	//mat1.print();
	Matrix mat2(size);
	mat2.set();
	//cout << mat1.det() << endl;
	Matrix mat3(size);
	int beg = clock();
	mat3 = multiStrassen(mat1, mat2, 0);
	int end = clock();
	cout << endl << endl << end - beg << endl;
	//beg = clock();
	////mat1 = mat1 * mat2;
	//end = clock();
	//cout << endl << endl << end - beg << endl;
	//mat3.print();
	system("PAUSE");
}