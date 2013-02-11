#ifndef _GOPT_LBFGS_H_
#define _GOPT_LBFGS_H_

#include <deque>
#include "array.h"

//
// LBFGS method
// following Nocedal
//
template <class LineSearch>
class LBFGSSolver
{
public:
    typedef typename LineSearch::Energy Energy;
    
protected:
    // interface
    bool m_debug;
    Energy* m_energy;
    
    // LBFGS internal data
    unsigned m_size;
    std::deque<double> m_rho;
    std::deque<Array> m_dx; // dx[k] =    x[k+1] -    x[k]
    std::deque<Array> m_dg; // dg[k] = grad[k+1] - grad[k]
    
public:
    LBFGSSolver(Energy* energy, bool debug, unsigned size)
    {
        m_energy = energy;
        m_debug = debug;
        m_size = size;
    }
    
    void clear()
    {
        m_dg.clear();
        m_dx.clear();
        m_rho.clear();
    }
    
    bool get_debug() const
    {
        return m_debug;
    }
    
    void set_debug(bool debug)
    {
        m_debug = debug;
    }
    
    unsigned get_size_limit() const
    {
        return m_size;
    }
    
    void set_size_limit(unsigned nb)
    {
        clear();
        m_size = nb;
    }
    
    unsigned get_current_size() const
    {
        return m_rho.size();
    }
    
    unsigned run(const unsigned max_iters,
                 const double tolerance,            // min gradient norm
                 const int  step_control,
                 const double step_size,
                 const unsigned max_search_iters) // # line searches
    {
        Array x0 = m_energy->get_variables();
        Array g0 = m_energy->gradient();
        
        unsigned iters = 0;
        while (iters < max_iters)
        {
            double norm = g0.length();
            if (m_debug) std::cout << "Norm: " << norm << std::endl;
            if (norm < tolerance) break;
            
            Array delta = solve_quasi_newton(g0);
            m_energy->project(delta);
            
            LineSearch search(m_energy, m_debug);
            bool ok = search.run(x0, delta, step_control, step_size, max_search_iters);
            if (!ok) break;
            iters++;
            
            Array x1 = m_energy->get_variables();
            Array g1 = m_energy->gradient();
            append(x1-x0, g1-g0);
            x0 = x1;
            g0 = g1;
        }
        return iters;
    }
    
    Array solve_quasi_newton(const Array& g)
    {
        Array delta = g;
        unsigned nb = get_current_size();
        if (nb > 0)
        {
            // reverse traverse
            std::vector<double> alpha(nb);
            for (unsigned i = 0; i < nb; ++i)
            {
                unsigned index = nb - 1 - i;
                alpha[index] = m_rho[index] * dot(delta, m_dx[index]);
                delta -= alpha[index] * m_dg[index];
            }
            
            delta *= compute_gamma();
            
            // forward traverse
            for (unsigned i = 0; i < nb; ++i)
            {
                unsigned index = i;
                double beta = m_rho[index] * dot(delta, m_dg[index]);
                delta += (alpha[index] - beta) * m_dx[index];
            }
        }
        return -delta;
    }
    
    void append(const Array& dx, const Array& dg)
    {
        if (get_size_limit() == 0) return;
        
        double rho = dot(dx, dg);
        if (rho != 0.0) rho = 1.0/rho;
        
        m_rho.push_back(rho);
        m_dx.push_back(dx);
        m_dg.push_back(dg);
        
        if (get_current_size() > get_size_limit())
        {
            m_rho.pop_front();
            m_dx.pop_front();
            m_dg.pop_front();
        }
    }
    
    double compute_gamma()
    {
        if (get_current_size() == 0) return 1.0;
        
        double D = m_dg.back().length2();
        if (D == 0.0) return 1.0;
        
        double N = dot(m_dx.back(), m_dg.back());
        return (N / D);
    }
};

#endif
