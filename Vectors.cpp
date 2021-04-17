#include "Vectors.h"


Vector::Vector(int dimension) : size(dimension)
{
	elements = new double[size];
	for (int i = 0; i < size; ++i)
	{
		elements[i] = 0.;
	}
}

Vector::Vector(const Vector& vector) : size(vector.size)
{
	elements = new double[size];
	for (int i = 0; i < size; ++i)
	{
		elements[i] = vector[i];
	}
}

void Vector::addElement(int position, double element)
{
	elements[position] += element;
}

double Vector::findNorm() const
{
	double result = 0;
	for (int i = 0; i < size; ++i)
	{
		result += elements[i] * elements[i];
	}
	return sqrt(result);
}

Vector& Vector::operator =(const Vector& vector)
{
	size = vector.getSize();
	for (int i = 0; i < size; ++i)
	{
		elements[i] = vector[i];
	}
	return (*this);
}

Vector Vector::operator *(const double& number) const
{
	Vector result((*this));
	for (int i = 0; i < size; ++i)
	{
		result[i] *= number;
	}
	return result;
}

double Vector::operator *(const Vector& vector) const
{
	double result = 0;
	for (int i = 0; i < size; ++i)
	{
		result += elements[i] * vector[i];
	}
	return result;
}

Vector& Vector::operator +=(const Vector& vector)
{
	for (int i = 0; i < size; ++i)
	{
		elements[i] += vector[i];
	}
	return (*this);
}

Vector& Vector::operator -=(const Vector& vector)
{
	for (int i = 0; i < size; ++i)
	{
		elements[i] -= vector[i];
	}
	return (*this);
}

Vector Vector::operator +(const Vector& vector) const
{
	Vector result((*this).getSize());
	for (int i = 0; i < size; ++i)
	{
		result[i] = elements[i] + vector[i];
	}
	return result;
}

Vector Vector::operator -(const Vector& vector) const
{
	Vector result((*this).getSize());
	for (int i = 0; i < size; ++i)
	{
		result[i] = elements[i] - vector[i];
	}
	return result;
}

std::ostream& operator <<(std::ostream& out_stream, const Vector& vector)
{
	for (int i = 0; i < vector.size; ++i)
	{
		out_stream << vector[i] << " " << "\n";
	}
	return out_stream;
}