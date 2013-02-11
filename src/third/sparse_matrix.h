#ifndef _SPARSE_MATRIX_H_
#define _SPARSE_MATRIX_H_

/*
 * class Sparsematrix
 * any rectangular sparse matrix
 * values stored by row (CRS)
 * CRS = Compressed by Row Storage
 */

#include "sparse_array.h"

class SparseMatrix {
private:
    std::vector<SparseArray> _row;
    
public:
    // Ctor
    SparseMatrix(unsigned nrows = 1, unsigned ncols = 1);
    
    SparseMatrix(const SparseMatrix& rhs);
    
    // Dtor
    ~SparseMatrix();
    
    // Access
    
    unsigned numRows() const;
    
    unsigned numColumns() const;
    
    unsigned numNonZeros() const;
    
    unsigned numLowerNonZeros() const;
    
    unsigned copyCRS(long*& rowptr, long*& index, double*& value) const;
    
    unsigned copyCCS(long*& colptr, long*& index, double*& value) const;
    
    double getValue(unsigned i, unsigned j) const;
    
    void setValue(unsigned i, unsigned j, const double value);
    
    SparseArray getRow(unsigned i) const;
    
    SparseArray getColumn(unsigned i) const;
    
    SparseMatrix getTranspose() const;
    
    void setRow(unsigned i, const SparseArray& row);
    
    void setColumn(unsigned i, const SparseArray& col);
    
    // Operators
    
    SparseMatrix& operator = (const SparseMatrix& rhs);
    
    SparseMatrix& operator += (const SparseMatrix& rhs);
    
    SparseMatrix& operator -= (const SparseMatrix& rhs);
    
    SparseMatrix& operator *= (const SparseMatrix& rhs);
    
    SparseMatrix& operator *= (const double rhs);
    
    SparseMatrix& operator /= (const double rhs);
    
    SparseMatrix multiply(const double b) const;
    
    SparseMatrix multiply(const SparseMatrix& b) const;
};

// Friend Operators

inline
SparseMatrix operator - (const SparseMatrix& rhs) {
    SparseMatrix A(rhs.numRows(), rhs.numColumns());
    for (unsigned i=0; i<rhs.numRows(); ++i) {
        A.setRow(i, -rhs.getRow(i));
    }
    return A;
}

inline
SparseMatrix operator + (const SparseMatrix& a, const SparseMatrix& b) {
    SparseMatrix c = a;
    c += b;
    return c;
}

inline
SparseMatrix operator - (const SparseMatrix& a, const SparseMatrix& b) {
    SparseMatrix c = a;
    c -= b;
    return c;
}

inline
SparseMatrix operator * (const SparseMatrix& a, const SparseMatrix& b) {
    return a.multiply(b);
}

inline
SparseMatrix operator * (const SparseMatrix& a, const double b) {
    return a.multiply(b);
}

inline
SparseMatrix operator * (const double a, const SparseMatrix& b) {
    return b.multiply(a);
}

inline
SparseMatrix operator / (const SparseMatrix& a, const double b) {
    SparseMatrix c(a.numRows(), a.numColumns());
    for (unsigned i=0; i<a.numRows(); ++i)
    {
        c.setRow(i, a.getRow(i) / b);
    }
    return c;
}

inline
Array operator * (const SparseMatrix& M, const Array& x) {
    Array y(M.numRows());
    for (unsigned i=0; i<M.numRows(); ++i)
    {
        y[i] = dot(M.getRow(i), x);
    }
    return y;
}

inline
Array operator * (const Array& x, const SparseMatrix& M) {
    Array y(M.numColumns());
    for (unsigned i = 0; i < M.numRows(); ++i)
    {
        const SparseArray& row = M.getRow(i);
        for (unsigned j = 0; j < row.numNonZeros(); ++j)
        {
            unsigned k = row.readIndex(j);
            double Mik = row.readValue(j);
            y[k] += x[i]*Mik;
        }
    }
    return y;
}

// Concatenation

inline
SparseMatrix ConcatenateV(const SparseMatrix& A, const SparseMatrix& B)
{
    SparseMatrix C(A.numRows() + B.numRows(), A.numColumns());
    
    // copy A
    for (unsigned i=0; i<A.numRows(); ++i) {
        C.setRow(i, A.getRow(i));
    }
    
    // copy B
    for (unsigned i=0; i<B.numRows(); ++i) {
        C.setRow(A.numRows() + i, B.getRow(i));
    }
    
    return C;
}

inline
SparseMatrix ConcatenateH(const SparseMatrix& A, const SparseMatrix& B)
{
    unsigned numColumns = A.numColumns() + B.numColumns();
    SparseMatrix C(A.numRows(), numColumns);
    
    for (unsigned i = 0; i < C.numRows(); ++i)
    {
        SparseArray Crow(numColumns);
        
        const SparseArray& Arow = A.getRow(i);
        for (unsigned k = 0; k < Arow.numNonZeros(); ++k)
        {
            unsigned index = Arow.readIndex(k);
            double value = Arow.readValue(k);
            Crow.setValue(index, value);
        }
        
        const SparseArray& Brow = B.getRow(i);
        for (unsigned k = 0; k < Brow.numNonZeros(); ++k)
        {
            unsigned index = A.numColumns() + Brow.readIndex(k);
            double value = Brow.readValue(k);
            Crow.setValue(index, value);
        }
        
        C.setRow(i, Crow);
    }
    
    return C;
}

#endif
