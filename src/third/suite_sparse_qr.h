#ifndef _SUITESPARSEQR_WRAPPER_H
#define _SUITESPARSEQR_WRAPPER_H 1

// SuiteSparseQR header
#include <SuiteSparseQR.hpp>
#include <cholmod.h>

// local
#include "array.h"
#include "sparse_matrix.h"

class SuiteSparseQRFactorizer
{
private:
    cholmod_common _common;
    SuiteSparseQR_factorization<double>* _QR;
    
public:
    SuiteSparseQRFactorizer()
    {
        cholmod_l_start(&_common);
        _QR = NULL;
    }
    
    ~SuiteSparseQRFactorizer()
    {
        if (_QR) SuiteSparseQR_free<double>(&_QR, &_common);
        cholmod_l_finish(&_common);
    }
    
    ////////////////////////////////////////////
    // solve (under)-determined linear system //
    ////////////////////////////////////////////
    
    bool factorize(const SparseMatrix& A)
    {
        if (_QR) return false;
        
        // A'*E = Q*R -> A = E*R'*Q', inv(E) = E'
        cholmod_sparse* At = convert_transpose_to_cholmod_sparse(A);
        _QR = SuiteSparseQR_factorize<double>(SPQR_ORDERING_DEFAULT,
                                              SPQR_DEFAULT_TOL,
                                              At,
                                              &_common);
        cholmod_l_free_sparse(&At, &_common);
        return true;
    }
    
    bool solve(const Array& rhs, Array& x)
    {
        if (!_QR) return false;
        
        // R'*y = E'*b
        // x = Q*y
        // A*x = b -> (E*R'*Q')*(Q*inv(R')*E'*b) = b
        
        cholmod_dense* b = convert_to_cholmod_dense(rhs);
        cholmod_dense* y = SuiteSparseQR_solve(SPQR_RTX_EQUALS_ETB, _QR, b, &_common);
        cholmod_dense* z = SuiteSparseQR_qmult(SPQR_QX, _QR, y, &_common);
        x = convert_to_array(z);
        
        cholmod_l_free_dense(&b, &_common);
        cholmod_l_free_dense(&y, &_common);
        cholmod_l_free_dense(&z, &_common);
        return true;
    }
    
    ///////////////////////////////////////////////
    // Solve under/over-determined linear system //
    ///////////////////////////////////////////////
    
    void solve_once(const SparseMatrix& A,
                    Array& x,
                    const Array& b)
    {
        cholmod_sparse* matrix = convert_to_cholmod_sparse(A);
        cholmod_dense* rhs = convert_to_cholmod_dense(b);
        cholmod_dense* solution = SuiteSparseQR<double>(matrix, rhs, &_common);
        x = convert_to_array(solution);
        cholmod_l_free_dense(&solution, &_common);
        cholmod_l_free_dense(&rhs, &_common);
        cholmod_l_free_sparse(&matrix, &_common);
    }
    
    /////////////////////////////////////////
    // solve over-determined linear system //
    /////////////////////////////////////////
    
    bool factorize_for_least_squares(const SparseMatrix& A)
    {
        if (_QR) return false;
        
        // A*E = Q*R -> A = Q*R*E', inv(E) = E'
        cholmod_sparse* matrix = convert_to_cholmod_sparse(A);
        _QR = SuiteSparseQR_factorize<double>(SPQR_ORDERING_DEFAULT,
                                              SPQR_DEFAULT_TOL,
                                              matrix,
                                              &_common);
        cholmod_l_free_sparse(&matrix, &_common);
        return true;
    }
    
    bool solve_for_least_squares(const Array& rhs, Array& x)
    {
        if (!_QR) return false;
        
        // y = Q'*b
        // R*E'*x = y
        // R*E'*x = Q'*b -> Q*R*E'*x = b -> A*x = b
        cholmod_dense* b = convert_to_cholmod_dense(rhs);
        cholmod_dense* y = SuiteSparseQR_qmult(SPQR_QTX, _QR, b, &_common);
        cholmod_dense* z = SuiteSparseQR_solve(SPQR_RETX_EQUALS_B, _QR, y, &_common);
        x = convert_to_array(z);
        
        cholmod_l_free_dense(&b, &_common);
        cholmod_l_free_dense(&y, &_common);
        cholmod_l_free_dense(&z, &_common);
        return true;
    }
    
private:
    cholmod_sparse* convert_to_cholmod_sparse(const SparseMatrix& A)
    {
        // B = A
        long* colptr;
        long* index;
        double* value;
        int nnz = (int) A.copyCCS(colptr, index, value);
        cholmod_sparse* B = cholmod_l_allocate_sparse(A.numRows(), A.numColumns(), nnz,
                                                      1, 1, 0, CHOLMOD_REAL, &_common);
        B->p = (void*) colptr;
        B->i = (void*) index;
        B->x = (void*) value;
        return B;
    }
    
    cholmod_sparse* convert_transpose_to_cholmod_sparse(const SparseMatrix& A)
    {
        // B = A'
        long* colptr;
        long* index;
        double* value;
        int nnz = (int) A.copyCRS(colptr, index, value);
        cholmod_sparse* B = cholmod_l_allocate_sparse(A.numColumns(), A.numRows(), nnz,
                                                      1, 1, 0, CHOLMOD_REAL, &_common);
        B->p = (void*) colptr;
        B->i = (void*) index;
        B->x = (void*) value;
        return B;
    }
    
    cholmod_dense* convert_to_cholmod_dense(const Array& a)
    {
        cholmod_dense* b = cholmod_l_allocate_dense(a.size(), 1, a.size(), CHOLMOD_REAL, &_common);
        double* tmp = (double*) b->x;
        for (unsigned i = 0; i < a.size(); ++i) tmp[i] = a[i];
        return b;
    }
    
    Array convert_to_array(cholmod_dense* a)
    {
        Array b(a->nrow);
        double* tmp = (double*) a->x;
        for (unsigned i = 0; i < a->nrow; ++i) b[i] = tmp[i];
        return b;
    }
};

#endif
