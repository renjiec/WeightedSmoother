#ifndef _GOPT_NEWTON_SOLVER_H_
#define _GOPT_NEWTON_SOLVER_H_

#include "array.h"
#include "sparse_matrix.h"

//
// Newton's method
//
template <class LineSearch, class LinearSolver>
class NewtonSolver
{
public:
    typedef typename LineSearch::Energy Energy;
    
protected:
    Energy* m_energy;
    bool m_debug;
    
public:
    NewtonSolver(Energy* energy, bool debug)
    {
        m_energy = energy;
        m_debug = debug;
    }
    
    unsigned run(const unsigned max_iters,        // # Newton steps
                 const double tolerance,          // min gradient norm
                 const int  step_control,
                 const double step_size,
                 const unsigned max_search_iters, // # line searches
                 const double perturbation)       // Hessian + pert*I
    {
        unsigned iters = 0;
        while (iters < max_iters)
        {
            Array grad = m_energy->gradient();
            double norm = grad.length();
            if (m_debug) std::cout << "Norm: " << norm << std::endl;
            if (norm < tolerance) break;
            
            SparseMatrix* H = m_energy->hessian();
            if (perturbation != 0.0) regularize_matrix(*H, perturbation*norm);
            
            Array delta = solve_linear_system(*H, -grad);
            m_energy->project(delta);
            delete(H);
            
            LineSearch search(m_energy, m_debug);
            Array x = m_energy->get_variables();
            bool ok = search.run(x, delta, step_control, step_size, max_search_iters);
            if (!ok) break;
            iters++;
        }
        return iters;
    }
    
    void regularize_matrix(SparseMatrix& A, const double epsilon) const
    {
        for (unsigned i = 0; i < A.numRows(); ++i)
        {
            double Aii = A.getValue(i, i);
            A.setValue(i, i, Aii + epsilon);
        }
    }
    
    Array solve_linear_system(const SparseMatrix& A,
                                     const Array& b) const
    {
        Array x;
        LinearSolver solver;
        solver.solve_once(A, x, b);
        return x;
    }
};

#endif
