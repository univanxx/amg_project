#ifndef VECTORS_H
#define VECTORS_H

#include <math.h>
#include <iostream>

class Vector
{
private:
	int size;
	double* elements;

public:
	Vector(int dimension);  // �������� ������� ������ �������� ����� 
	Vector(const Vector& vector); // �����������
	int getSize() const  // ������ �������
	{
		return size;
	};
	double* getElements() const  // ��������� ��������� �������
	{
		return elements;
	};
	void addElement(int position, double element);  // ��������� ������ � �������� � ������ ��������
	double findNorm() const;  // ����� L_{2}
	Vector& operator =(const Vector& vector);  // �������� ������������
	Vector& operator =(const double* vector);  // ��� ���� ������� ������������
	Vector operator *(const double& number) const;  // �������� ��������� �� ������
	Vector operator /(const double& number) const;  // �������� ������� �� ������
	double operator *(const Vector& vector) const;  // �������� ��������� �� ������
	Vector& operator +=(const  Vector& vector);  // �������� �������� ��������
	Vector& operator -=(const  Vector& vector);  // �������� ��������� ��������
	Vector operator +(const  Vector& vector) const;  // ��� ���� �������� �������� ��������
	Vector operator -(const  Vector& vector) const;  // ��� ���� �������� ���������� ��������
	double& operator [](int index) const  // �������� ������ �������
	{
		return elements[index];
	};
	friend std::ostream& operator<<(std::ostream&, const Vector&); // ������� ������ �������
	~Vector()  // ����������
	{
		delete[] elements;
	};
};


#endif

