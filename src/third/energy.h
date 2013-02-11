#ifndef _TEMPLATE_ENERGY_H_
#define _TEMPLATE_ENERGY_H_

#include "array.h"
#include "sparse_matrix.h"

class TemplateEnergy
{
public:
    virtual ~TemplateEnergy() { }
    
    virtual Array get_variables() const = 0;
    
    virtual bool set_variables(const Array& x) = 0;
    
    virtual double evaluate() const = 0;
    
    virtual Array gradient() const = 0;
    
    virtual SparseMatrix* hessian() const { return NULL; }
    
    virtual void project(Array&) const { }
};

#endif
