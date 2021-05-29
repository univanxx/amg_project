//#include "stdafx.h"
#include "Matrices.h"
#include "Vectors.h"
Matrix::Matrix(int number_rows, int number_columns) : rows(number_rows), columns(number_columns)
{
	elements = new double* [rows];
	for (int i = 0; i < rows; ++i)
	{
		elements[i] = new double[columns];
		for (int j = 0; j < columns; ++j)
		{
			elements[i][j] = 0.;
		}
	}
}

Matrix::Matrix(const Matrix& matrix) : rows(matrix.rows), columns(matrix.columns)
{
	elements = new double* [rows];
	for (int i = 0; i < rows; ++i)
	{
		elements[i] = new double[columns];
		for (int j = 0; j < columns; ++j)
		{
			elements[i][j] = matrix.getElements()[i][j];
		}
	}
}

void Matrix::setElement(int row, int column, double element)
{
	(*this).elements[row][column] = element;
}

void Matrix::addElement(int row, int column, double element)
{
	(*this).elements[row][column] += element;
}

void Matrix::addVector(int row, const Vector& vector)
{
	for (int i = 0; i < (*this).getColumns(); ++i)
	{
		(*this).addElement(row, i, vector.getElements()[i]);
	}
}
Matrix Matrix::transpose() const
{
	Matrix result(columns, rows);
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < columns; ++j)
		{
			result.setElement(j, i, (*this).getElements()[i][j]);
		}
	}
	return result;
}

Matrix& Matrix::operator =(const Matrix& matrix)
{
	rows = matrix.getRows();
	columns = matrix.getColumns();
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < columns; ++j)
		{
			elements[i][j] = matrix.getElements()[i][j];
		}
	}
	return (*this);
}

Matrix Matrix::operator *(const Matrix& matrix)
{
	Matrix result(rows, matrix.columns);
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < matrix.columns; ++j)
		{
			for (int k = 0; k < columns; ++k)
			{
				result.addElement(i, j, (*this).getElements()[i][k] * matrix.getElements()[k][j]);
			}
		}

	}
	return result;
}

Vector Matrix::operator [](const int index)
{
	int len = (*this).getColumns();
	Vector res(len);
	for (int k = 0; k < len; ++k)
	{
		res[k] = (*this).getElements()[index][k];
	}
	return res;
}

Vector Matrix::operator *(const Vector& vector) const
{
	Vector result(rows);
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < columns; ++j)
		{
			result[i] += (*this).getElements()[i][j] * vector[j];
		}
	}
	return result;
}

Matrix Matrix::operator *(const double& number)
{
	Matrix result(*this);
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < columns; ++j)
		{
			result.setElement(i, j, (*this).getElements()[i][j] * number);
		}
	}
	return  result;
}

Matrix Matrix::operator +(const Matrix& matrix)
{
	Matrix result(rows, matrix.columns);
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < matrix.columns; ++j)
		{
			result.setElement(i, j, (*this).getElements()[i][j] + matrix.getElements()[i][j]);
		}
	}
	return result;
}

std::ostream& operator<<(std::ostream& out_stream, const Matrix& matrix)
{
	for (int i = 0; i < matrix.rows; ++i)
	{
		for (int j = 0; j < matrix.columns; ++j)
		{
			out_stream << matrix.getElements()[i][j] << "  ";
		}
		out_stream << "\n";
	}
	return out_stream;
}