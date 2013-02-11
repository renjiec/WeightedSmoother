#ifndef _PHALFEDGE_H_
#define _PHALFEDGE_H_

#include "dpqueue.h"

template <class FT, class Halfedge_handle>
class CPHalfedge
{
protected:
    Halfedge_handle m_halfedge;
    FT              m_priority;
    
public:
    CPHalfedge()
    {
        m_halfedge = Halfedge_handle();
        m_priority = 0.0;
    }
    
    CPHalfedge(const Halfedge_handle& he, const FT priority = 0.0)
    {
        m_halfedge = he;
        m_priority = priority;
    }
    
    CPHalfedge(const CPHalfedge& phedge) 
    {
        m_halfedge = phedge.halfedge();
        m_priority = phedge.priority();
    }
    
    virtual ~CPHalfedge() { }
    
    CPHalfedge& operator = (const CPHalfedge& phedge) 
    {
        m_halfedge = phedge.halfedge();
        m_priority = phedge.priority();
        return *this;
    }
    
    bool operator == (const CPHalfedge& phedge) const
    {
        return m_halfedge == phedge.halfedge();
    }
    
    bool operator < (const CPHalfedge& phedge) const
    {
			return &*m_halfedge < &*phedge.halfedge();
    }
    
    const Halfedge_handle halfedge() const { return m_halfedge; }
    const FT priority() const { return m_priority; }
};

template <class T>
class CDPQueue_short : public DynamicPriorityQueue<T>
{
public:
    CDPQueue_short() { }
    
    ~CDPQueue_short() { }
    
    bool compare(const T& a, const T& b) const
    {
        return a.priority() < b.priority();
    }
};

template <class T>
class CDPQueue_long : public DynamicPriorityQueue<T>
{
public:
    CDPQueue_long() { }
    
    ~CDPQueue_long() { }
    
    bool compare(const T& a, const T& b) const
    {
        return a.priority() > b.priority();
    }
};

#endif
