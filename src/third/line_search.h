#ifndef _GOPT_LINE_SEARCH_H_
#define _GOPT_LINE_SEARCH_H_

/*
 * Based on "Numerical Optimizaiton"
 * by J. Nocedal and S. Wright
 * pages 60-61
 *
 * phi(t) = E[x0 + t*v]
 * D[phi,t] = D[E,x]*v
 *
 */

#include "array.h"

template <class CENERGY>
class LineSearch
{
public:
    typedef CENERGY Energy;
    
    enum
    {
        NONE = 0,
        ENERGY_BASED = 1,
        GRAD_BASED = 2,
        VALIDITY_BASED = 3
    };
    
protected:
    Energy* m_energy;
    bool m_debug;
    
public:
    LineSearch(Energy* energy, bool debug = false)
    {
        m_energy = energy;
        m_debug = debug;
    }
    
    bool run(const Array& x0,
             const Array& v,
             const int control,
             const double step,
             const unsigned max_iters)
    {
        if (control == NONE)
            return move(x0, v, step);
        
        if (control == ENERGY_BASED)
            return energy_back_track(x0, v, step, max_iters);
        
        if (control == GRAD_BASED)
            return grad_back_track(x0, v, step, max_iters);
        
        if (control == VALIDITY_BASED)
            return validity_back_track(x0, v, step, max_iters);
        
        return false;
    }
    
    bool move(const Array& x0,
              const Array& v,
              const double step)
    {
        bool ok = try_move(x0, v, step);
        if (m_debug) std::cout << "LS:Move" << std::endl;
        if (ok) return true;
        
        // reset
        if (m_debug) std::cout << "LS:Reset" << std::endl;
        m_energy->set_variables(x0);
        return false;
    }
    
    bool validity_back_track(const Array& x0,
                             const Array& v,
                             const double max_step,
                             const unsigned max_iters)
    {
        double lower_alpha = 0.0;
        double upper_alpha = max_step;
        
        for (unsigned i = 0; i < max_iters; ++i)
        {
            double alpha = pick_alpha(lower_alpha, upper_alpha);
            bool ok = try_move(x0, v, alpha);
            if (ok) return true;
            upper_alpha = alpha;
        }
        
        // reset
        if (m_debug) std::cout << "LS:Reset" << std::endl;
        m_energy->set_variables(x0);
        return false;
    }
    
    bool energy_back_track(const Array& x0,
                           const Array& v,
                           const double max_step,
                           const unsigned max_iters,
                           const double c1 = 1.0e-4)
    {
        double phi0  = m_energy->evaluate();
        double gphi0 = dot(v, m_energy->gradient());
        
        double lower_alpha = 0.0;
        double upper_alpha = max_step;
        
        for (unsigned i = 0; i < max_iters; ++i)
        {
            double alpha = pick_alpha(lower_alpha, upper_alpha);
            
            bool ok = try_move(x0, v, alpha);
            if (ok)
            {
                double phi = m_energy->evaluate();
                double linear_approx = phi0 + c1 * alpha * gphi0;
                
                /*
                 if (m_debug)
                 {
                 std::cout << "LS::Phi: " << phi << std::endl;
                 std::cout << "LS::LinearApprox: " << linear_approx << std::endl;
                 }
                 */
                
                if (phi <= linear_approx)
                {
                    if (m_debug) std::cout << "LS::FinalAlpha = " << alpha << std::endl;
                    return true;
                }
            }
            upper_alpha = alpha;
        }
        
        // reset
        if (m_debug) std::cout << "LS:Reset" << std::endl;
        m_energy->set_variables(x0);
        return false;
    }
    
    bool grad_back_track(const Array& x0,
                         const Array& v,
                         const double max_step,
                         const unsigned max_iters)
    {
        Array g0 = m_energy->gradient();
        double max_g0 = std::max( abs(g0.getMin()), abs(g0.getMax()) );
        double norm0 = g0.length();
        
        double lower_alpha = 0.0;
        double upper_alpha = max_step;
        
        for (unsigned i = 0; i < max_iters; ++i)
        {
            double alpha = pick_alpha(lower_alpha, upper_alpha);
            
            bool ok = try_move(x0, v, alpha);
            if (ok)
            {
                Array g = m_energy->gradient();
                double max_g = std::max( abs(g.getMin()), abs(g.getMax()) );
                double norm = g.length();
                
                //if (max_g < max_g0)
                if (norm <= norm0)
                {
                    if (m_debug) std::cout << "LS::FinalAlpha = " << alpha << std::endl;
                    return true;
                }
            }
            upper_alpha = alpha;
        }
        
        // reset
        if (m_debug) std::cout << "LS:Reset" << std::endl;
        m_energy->set_variables(x0);
        return false;
    }
    
protected:
    double pick_alpha(const double lower, const double upper) const
    {
        return 0.5*(lower + upper);
    }
    
    bool try_move(const Array& x0,
                  const Array& v,
                  const double step)
    {
        Array x = x0 + step*v;
        bool ok = m_energy->set_variables(x);
        if (m_debug) std::cout << "LS::try: step = " << step << " -> ok = " << ok << std::endl;
        return ok;
    }
};

#endif
