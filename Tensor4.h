#ifndef TENSOR4_H
#define TENSOR4_H

#include "Tensors.h"
#include <iostream>

class Tensor4
{
private:
	int axle1;
	int axle2;
	int rows;
	int columns;
	double**** elements;

public:
	Tensor4(int number_axis1, int number_axis2, int number_rows, int number_columns);  // �������� ������� ������ �������� ��������
	Tensor4(const Tensor4& tensor);  // �����������
	double**** getElements() const  // ��������� ��������� �������
	{
		return elements;
	};
	int getAxle1() const  // ��������� ���������� ����
	{
		return axle1;
	};
	int getAxle2() const  // ��������� ���������� ����
	{
		return axle2;
	};
	int getRows() const  // ��������� ���������� �����
	{
		return rows;
	};
	int getColumns() const  // ��������� ���������� ��������
	{
		return columns;
	};
	void setElement(int axle1, int axle2, int row, int column, double element);  // �������� �������� � �������� �������
	void addElement(int axle1, int axle2, int row, int column, double element);  // ���������� �������� � ��������� ��������
	void setMatrix(int axle1, int axle2, const Matrix &M); // ���������� ������� �� �������
	void addMatrix(int axle1, int axle2, const Matrix &M); // ���������� ������� �� �������

	Tensor4& operator =(const Tensor4& tensor);  // �������� ������������
	Tensor3 operator [](const int index);

	friend std::ostream& operator<<(std::ostream&, const Tensor4& tensor);
	~Tensor4()  // ����������
	{
		for (int m = 0; m < axle1; ++m)
		{
			for (int i = 0; i < axle2; ++i)
			{
				for (int j = 0; j < rows; ++j)
				{
					delete[] elements[m][i][j];
				}
				delete elements[m][i];
			}
			delete[] elements[m];
		}
		delete[] elements;
	};
};


#endif
