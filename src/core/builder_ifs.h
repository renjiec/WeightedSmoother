#ifndef _BUILDER_IFS_
#define _BUILDER_IFS_

#include <CGAL/Polyhedron_incremental_builder_3.h>

template <
class HDS,
class Kernel,
class Point_iterator,
class Index_iterator
>
class CModifierIfs : public CGAL::Modifier_base<HDS>
{
public:
    typedef typename Kernel::Point_2      Point_2;
    typedef typename Kernel::Point_3      Point_3;
    typedef typename HDS::Face_handle     Face_handle;
    typedef typename HDS::Halfedge_handle Halfedge_handle;
    typedef typename CGAL::Polyhedron_incremental_builder_3<HDS> Builder;
    
private:
	Point_iterator m_pbegin, m_pend;
	Index_iterator m_ibegin, m_iend;
    
public:
    CModifierIfs(Point_iterator pbegin, Point_iterator pend,
	             Index_iterator ibegin, Index_iterator iend)
    {
		m_pend = pend;
		m_iend = iend;
		m_pbegin = pbegin;
		m_ibegin = ibegin;
    }

    ~CModifierIfs() { }
    
    void operator()(HDS& hds)
    {
        Builder B(hds,true);
        B.begin_surface(3,1,6);
        add_vertices(B);
        add_facets(B);
        B.end_surface();
    }
    
    void add_vertices(Builder &B)
    {
        for (Point_iterator it = m_pbegin; it != m_pend; it++)
		{
			const Point_2& p = *it;
            B.add_vertex(Point_3(p.x(), p.y(), 0.0));
		}
    }
    
    void add_facets(Builder &B)
    {
	    for (Index_iterator it = m_ibegin; it != m_iend; )
        {
            B.begin_facet();
            B.add_vertex_to_facet(*it++);
            B.add_vertex_to_facet(*it++);
            B.add_vertex_to_facet(*it++);
            B.end_facet();
        }
    }
};

template <class Mesh>
class CBuilder_ifs
{
public:
    typedef typename Mesh::Kernel Kernel;
    typedef typename Mesh::HalfedgeDS HalfedgeDS;

    CBuilder_ifs() { }
    
    ~CBuilder_ifs() { }
    
    template <class Point_iterator, class Index_iterator>
    void run(Point_iterator pbegin, Point_iterator pend,
             Index_iterator ibegin, Index_iterator iend,
             Mesh &new_mesh)
    {
        CModifierIfs<HalfedgeDS, Kernel, Point_iterator, Index_iterator>
        builder(pbegin, pend, ibegin, iend);
        new_mesh.delegate(builder);
    }
};

#endif
