#ifndef _MESH_
#define _MESH_

#include <QtOpenGL>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/predicates/Regular_triangulation_ftC2.h>
#include <CGAL/convex_hull_2.h>

#include <vector>
#include <queue>

#include "random.h"
#include "primitives.h"
#include "../third/array.h"
#include "transparent.h"

#undef min
#undef max

// used by Dijkstra 
template<class Element>
struct Bigger_distance
{ 
  bool operator()(const Element& e1, 
                  const Element& e2) const
  {
    return e1->distance() > e2->distance();
  }
};


// compute facet normal
struct Facet_normal	// (functor)
{
	template <class	Facet>
	void operator()(Facet& f)
    {
        typename Facet::Normal sum = CGAL::NULL_VECTOR;
		typename Facet::Halfedge_around_facet_circulator h = f.facet_begin();
		do
		{
			typename Facet::Normal normal = CGAL::cross_product(h->next()->vertex()->point() - h->vertex()->point(),
                                                                h->next()->next()->vertex()->point() - h->next()->vertex()->point());
			double sqnorm = normal * normal;
            if (sqnorm != 0)
            {
                normal = normal / (double)std::sqrt(sqnorm);
            }
            sum	=	sum	+	normal;
        }
		while(++h	!= f.facet_begin());
        
		double	sqnorm = sum * sum;
		if(sqnorm	!= 0.0)
        {
			f.normal() = sum / std::sqrt(sqnorm);
        }
        else
        {
            f.normal() = CGAL::NULL_VECTOR;
            std::cout << "degenerated facet" << std::endl;
        }
    }
};
            
// compute vertex	normal
struct Vertex_normal //	(functor)
{
    template <class	Vertex>
    void operator()(Vertex&	v)
    {
        typename Vertex::Normal	normal = CGAL::NULL_VECTOR;
        typename Vertex::Halfedge_around_vertex_const_circulator	he	=	v.vertex_begin();
        typename Vertex::Halfedge_around_vertex_const_circulator	begin	=	he;
        CGAL_For_all(he,begin)
        if(!he->is_border())
            normal = normal	+	he->facet()->normal();
        double	sqnorm = normal * normal;
        if(sqnorm != 0.0)
            v.normal() = normal	/	(double)std::sqrt(sqnorm);
        else
            v.normal() = CGAL::NULL_VECTOR;
    }
};

template <class KERNEL, class ITEMS>
class CMesh : public CGAL::Polyhedron_3<KERNEL, ITEMS>
{
public:
	typedef KERNEL Kernel;
	typedef ITEMS  Items;
    
	typedef CGAL::Polyhedron_3<Kernel, Items> Base;
	typedef CMesh<Kernel, Items> Mesh;
	typedef typename Kernel::FT FT;
	// 2d
	typedef typename Kernel::Point_2    Point_2;
	typedef typename Kernel::Vector_2   Vector_2;
	typedef typename Kernel::Segment_2  Segment_2;
	typedef typename Kernel::Triangle_2 Triangle_2;
	typedef typename Kernel::Line_2     Line_2;
	typedef typename Kernel::Ray_2      Ray_2;
	// 3d
	typedef typename Kernel::Point_3    Point_3;
	typedef typename Kernel::Vector_3   Vector_3;
	typedef typename Kernel::Segment_3  Segment_3;
	typedef typename Kernel::Triangle_3 Triangle_3;
	typedef typename Kernel::Ray_3      Ray_3;
    
	typedef typename Mesh::Edge_iterator Edge_iterator;
	typedef typename Mesh::Facet_handle Facet_handle;
	typedef typename Mesh::Facet_iterator Facet_iterator;
	typedef typename Mesh::Vertex_handle Vertex_handle;
	typedef typename Mesh::Vertex_iterator Vertex_iterator;
	typedef typename Mesh::Halfedge_handle Halfedge_handle;
	typedef typename Mesh::Halfedge_iterator Halfedge_iterator;
	typedef typename Mesh::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;
	typedef typename Mesh::Halfedge_around_vertex_circulator Halfedge_around_vertex_circulator;
    
	enum { FIXED = 0, FREE = 1 };

	//-----------//
	// attributes //
	//-----------//
    
	// 2D convex hull
	std::list<Point_2> m_mesh_2d_ch;
    
    // harmonic coefficients
    std::vector<FT> m_alpha;
    
	// store generators to render them
	typedef std::vector<Halfedge_handle> Generator;
	std::vector<Generator> m_generators;
    
public:
	CMesh() { }
    
	~CMesh() { }

	std::vector<Generator>& generators() { return m_generators; }
    
	Vertex_handle get_source_vertex(Halfedge_handle h) const
	{
		return h->opposite()->vertex();
	}
    
	Vertex_handle get_target_vertex(Halfedge_handle h) const
	{
		return h->vertex();
	}
    
	Vertex_handle get_opposite_vertex(Halfedge_handle h) const
	{
		return h->next()->vertex();
	}
    
	bool is_vertex_grounded(Vertex_handle v) const
	{
		return (v->type() == FIXED);
	}

	void set_all_vertices_free()
	{
		for (Vertex_iterator
             v = this->vertices_begin();
             v != this->vertices_end();
             v++)
			v->type() = FREE;
	}

	void toggle_vertex_type(Vertex_handle v) 
	{
		if(v->type() == FIXED) 
			v->type() = FREE;
		else
			v->type() = FIXED;
	}
    
	bool is_edge_grounded(Halfedge_handle he) const
	{
		bool source = is_vertex_grounded(get_source_vertex(he));
		bool target = is_vertex_grounded(get_target_vertex(he));
		return (source && target);
	}
    
    bool is_free_boundary(Halfedge_handle he) const
    {
        if (!is_border_edge(he)) return false;
        if (is_edge_grounded(he)) return false;
        return true;
    }
    
	// CENTROID //
    
    Point_2 voronoi_centroid(Vertex_handle v, bool& boundary) const
	{
		if (is_border_vertex(v))
		{
			boundary = true;
			return boundary_barycenter_2d(v);
		}
		boundary = false;
		return inner_voronoi_centroid(v);
	}
    
	Point_3 barycenter_3d(Vertex_handle v, bool& boundary) const
	{
		if (is_border_vertex(v))
		{
			boundary = true;
			return boundary_barycenter_3d(v);
		}
		boundary = false;
		return inner_barycenter_3d(v);
	}
    
	Point_2 barycenter_2d(Vertex_handle v, bool& boundary) const
	{
		if (is_border_vertex(v))
		{
			boundary = true;
			return boundary_barycenter_2d(v);
		}
		boundary = false;
		return inner_barycenter_2d(v);
	}
    
	FT min_2d_incident_edge_length(Vertex_handle v)
	{
		FT min_len = 1e6;
		const Point_2 pivot = v->get_2d_point();
		Halfedge_around_vertex_circulator he = v->vertex_begin();
		Halfedge_around_vertex_circulator end = he;
		CGAL_For_all(he, end)
		{
			const Point_2 p = get_source_vertex(he)->get_2d_point();
			FT len = std::sqrt(CGAL::squared_distance(pivot, p));
			if (len < min_len) min_len = len;
		}
		return min_len;
	}
    
	FT average_2d_incident_edge_length(Vertex_handle v)
	{
		FT sum = 0.0;
		const Point_2 pivot = v->get_2d_point();
		Halfedge_around_vertex_circulator he = v->vertex_begin();
		Halfedge_around_vertex_circulator end = he;
		int degree = 0;
		CGAL_For_all(he, end)
		{
			const Point_2 p = get_source_vertex(he)->get_2d_point();
			FT len = std::sqrt(CGAL::squared_distance(pivot, p));
			sum += len;
			degree++;
		}
		return sum / (FT)degree;
	}
    
	Point_2 boundary_barycenter_2d(Vertex_handle v) const
	{
		int nb_neighbors = 0;
		const Point_2 pivot = v->get_2d_point();
		Halfedge_around_vertex_circulator he = v->vertex_begin();
		Halfedge_around_vertex_circulator end = he;
		Vector_2 vec = CGAL::NULL_VECTOR;
		CGAL_For_all(he, end)
		{
			// skip inner vertices
			// note that both he and opposite are required for checking
			if (!he->is_border() && !he->opposite()->is_border()) continue;
			const Point_2 neighbor = get_source_vertex(he)->get_2d_point();
			vec = vec + (neighbor - pivot);
			nb_neighbors++;
		}
		assert(nb_neighbors == 2);
		return pivot + (vec / double(nb_neighbors));
	}
    
	Point_2 inner_barycenter_2d(Vertex_handle v) const
	{
		int nb_neighbors = 0;
		const Point_2 pivot = v->get_2d_point();
		Halfedge_around_vertex_circulator vcir = v->vertex_begin();
		Halfedge_around_vertex_circulator end = vcir;
		Vector_2 vec = CGAL::NULL_VECTOR;
		CGAL_For_all(vcir, end)
		{
			const Point_2 neighbor = get_source_vertex(vcir)->get_2d_point();
			vec = vec + (neighbor - pivot);
			nb_neighbors++;
		}
		return pivot + (vec / double(nb_neighbors));
	}

	
	Point_3 odt_cc_3d(Vertex_handle v) const
	{
		Halfedge_around_vertex_circulator he = v->vertex_begin();
		Halfedge_around_vertex_circulator end = he;
		Vector_3 vec = CGAL::NULL_VECTOR;
		FT sum_areas = 0.0;
		CGAL_For_all(he, end)
		{
			if(he->is_border()) continue;

			Facet_handle f = he->facet();
			const Point_3 cc = circumcenter_3d(f);
			const FT facet_area = area_3d(f);
			vec = vec + facet_area * (cc - CGAL::ORIGIN);
			sum_areas += facet_area;
		}
		return CGAL::ORIGIN + vec / sum_areas;
	}

	Point_2 odt_cc_2d(Vertex_handle v) const
	{
		Halfedge_around_vertex_circulator he = v->vertex_begin();
		Halfedge_around_vertex_circulator end = he;
		Vector_2 vec = CGAL::NULL_VECTOR;
		FT sum_areas = 0.0;
		CGAL_For_all(he, end)
		{
			if(he->is_border()) continue;

			Facet_handle f = he->facet();
			const Point_2 cc = circumcenter(f);
			const FT facet_area = area(f);
			vec = vec + facet_area * (cc - CGAL::ORIGIN);
			sum_areas += facet_area;
		}
		return CGAL::ORIGIN + vec / sum_areas;
	}
    
	Point_3 boundary_barycenter_3d(Vertex_handle v) const
	{
		int nb_neighbors = 0;
		const Point_3 pivot = v->point();
		Halfedge_around_vertex_circulator he = v->vertex_begin();
		Halfedge_around_vertex_circulator end = he;
		Vector_3 vec = CGAL::NULL_VECTOR;
		CGAL_For_all(he, end)
		{
			// skip inner vertices
			// note that both he and opposite are required for checking
			if (!he->is_border() && !he->opposite()->is_border()) continue;
			const Point_3 neighbor = get_source_vertex(he)->point();
			vec = vec + (neighbor - pivot);
			nb_neighbors++;
		}
		assert(nb_neighbors == 2);
		return pivot + (vec / double(nb_neighbors));
	}
    
	Point_3 inner_barycenter_3d(Vertex_handle v) const
	{
		int nb_neighbors = 0;
		const Point_3 pivot = v->point();
		Halfedge_around_vertex_circulator vcir = v->vertex_begin();
		Halfedge_around_vertex_circulator end = vcir;
		Vector_3 vec = CGAL::NULL_VECTOR;
		CGAL_For_all(vcir, end)
		{
			const Point_3 neighbor = get_source_vertex(vcir)->point();
			vec = vec + (neighbor - pivot);
			nb_neighbors++;
		}
		return pivot + (vec / double(nb_neighbors));
	}
    
	Point_2 inner_voronoi_centroid(Vertex_handle v) const
	{
		// pick v as pivot to triangulate Voronoi cell
		const Point_2 pivot = v->get_2d_point();
		Vector_2 vec = CGAL::NULL_VECTOR;
		FT sum_areas = 0.0;
        
		// circulate around v
		Halfedge_around_vertex_circulator he = v->vertex_begin();
		Halfedge_around_vertex_circulator end = he;
		CGAL_For_all(he, end)
		{
			// assemble triangle (pivot c1 c2)
			Point_2 c1 = circumcenter(he->facet());
			Point_2 c2 = circumcenter(he->opposite()->facet());
			Triangle_2 triangle(pivot, c1, c2);
			const FT area = triangle.area();
			vec = vec + area * (CGAL::centroid(triangle) - CGAL::ORIGIN);
			sum_areas += area;
		}
		return CGAL::ORIGIN + (vec / sum_areas);
	}
    
	// BASIC GEOMETRY //
    
	FT area(Facet_handle f) const
	{
		return triangle(f).area();
	}
    
    FT area_3d(Facet_handle f) const
    {
        Triangle_3 t = triangle_3d(f);
        return std::sqrt(t.squared_area());
    }
    
    Triangle_2 triangle(Facet_handle f) const
	{
		Halfedge_handle he = f->halfedge();
		const Point_2 a = he->vertex()->get_2d_point();
		const Point_2 b = he->next()->vertex()->get_2d_point();
		const Point_2 c = he->next()->next()->vertex()->get_2d_point();
		return Triangle_2(a,b,c);
	}
    
    Triangle_3 triangle_3d(Facet_handle f) const
	{
		Halfedge_handle he = f->halfedge();
		const Point_3 a = he->vertex()->get_3d_point();
		const Point_3 b = he->next()->vertex()->get_3d_point();
		const Point_3 c = he->next()->next()->vertex()->get_3d_point();
		return Triangle_3(a,b,c);
	}
    
    Point_2 circumcenter(Facet_handle f) const
	{
		return CGAL::circumcenter(triangle(f));
	}

    Point_3 circumcenter_3d(Facet_handle f) const
	{
		return CGAL::circumcenter(triangle_3d(f));
	}
    
	Point_2 midpoint_2d(Halfedge_handle hij) const
	{
		const Point_2 pi = get_source_vertex(hij)->get_2d_point();
		const Point_2 pj = get_target_vertex(hij)->get_2d_point();
		return CGAL::midpoint(pi, pj);
	}
    
    Point_3 midpoint_3d(Halfedge_handle hij) const
	{
		const Point_3 pi = get_source_vertex(hij)->point();
		const Point_3 pj = get_target_vertex(hij)->point();
		return CGAL::midpoint(pi, pj);
	}
    
    Segment_2 segment_2d(Halfedge_handle hij) const
	{
		const Point_2 pi = get_source_vertex(hij)->get_2d_point();
		const Point_2 pj = get_target_vertex(hij)->get_2d_point();
		return Segment_2(pi, pj);
	}
    
	Segment_3 segment_3d(Halfedge_handle hij) const
	{
		const Point_3 pi = get_source_vertex(hij)->point();
		const Point_3 pj = get_target_vertex(hij)->point();
		return Segment_3(pi, pj);
	}
    
	FT squared_len_2d(Halfedge_handle hij) const
	{
		Segment_2 seg = segment_2d(hij);
		return seg.squared_length();
	}
    
	FT squared_len_3d(Halfedge_handle hij) const
	{
		Segment_3 seg = segment_3d(hij);
		return seg.squared_length();
	}
    
	FT len_2d(Halfedge_handle he) const
	{
		return std::sqrt(squared_len_2d(he));
	}
    
	FT len_3d(Halfedge_handle he) const
	{
		return std::sqrt(squared_len_3d(he));
	}
        
	Vector_2 get_vector_2d(Halfedge_handle hij) const
	{
		Segment_2 seg = segment_2d(hij);
		return seg.to_vector();
	}
    
    Vector_3 get_3d_vector(Halfedge_handle hij) const
	{
		Segment_3 seg = segment_3d(hij);
		return seg.to_vector();
	}
    
    Vector_2 rotate90(const Vector_2& vec) const
	{
		return Vector_2(-vec.y(), vec.x());
	}
    
	Vector_2 get_edge_normal(Halfedge_handle hij) const
	{
		Vector_2 eij = get_vector_2d(hij);
		FT lij = std::sqrt(eij*eij);
		if (lij != 0.0) eij = eij / lij;
		return rotate90(eij);
	}
    
	Vector_3 get_3d_edge_normal(Halfedge_handle hij) const
	{
        if (hij->is_border()) return Vector_3(0.0,0.0,0.0);
        Vector_3 n = get_3d_normal(hij->facet());
		Vector_3 eij = get_3d_vector(hij);
		FT lij = std::sqrt(eij*eij);
		if (lij != 0.0) eij = eij / lij;
		return CGAL::cross_product(n, eij);
	}
    
	Vector_3 get_3d_normal(Facet_handle f) const
	{
		Triangle_3 triangle = triangle_3d(f);
		Vector_3 n = triangle.supporting_plane().orthogonal_vector();
		FT norm = std::sqrt(n*n);
		if (norm != 0.0) n = n / norm;
		return n;
	}
    
	FT total_2d_area() const
	{
		FT sum_areas = 0.0;
		for (Facet_iterator
             f = this->facets_begin();
             f != this->facets_end();
             f++)
		{
			sum_areas += area(f);
		}
		return sum_areas;
	}
    
    FT total_3d_area()
	{
		FT sum_areas = 0.0;
		for (Facet_iterator
             f = this->facets_begin();
             f != this->facets_end();
             f++)
		{
			sum_areas += area_3d(f);
		}
		return sum_areas;
	}
    
    FT compute_min_len_3d()
	{
		FT min_len = 1e100;
		for (Edge_iterator
             he  = this->edges_begin();
             he != this->edges_end();
             he++)
		{
			min_len = std::min(min_len, len_3d(he));
		}
		return min_len;
	}
    
	FT compute_min_len_2d()
	{
		FT min_len = 1e100;
		for (Edge_iterator
             he  = this->edges_begin();
             he != this->edges_end();
             he++)
		{
			min_len = std::min(min_len, len_2d(he));
		}
		return min_len;
	}
    
	FT angle_2d(const Point_2& a, const Point_2& o, const Point_2& b) const
	{
		Vector_2 oa = a - o;
		Vector_2 ob = b - o;
		const FT norm_oa = std::sqrt(oa * oa);
		const FT norm_ob = std::sqrt(ob * ob);
		assert(norm_oa != 0.0);
		assert(norm_ob != 0.0);
		return acos((oa * ob) / (norm_oa * norm_ob));
	}
    
	FT angle_3d(const Point_3& a, const Point_3& o, const Point_3& b) const
	{
		Vector_3 oa = a - o;
		Vector_3 ob = b - o;
		const FT norm_oa = std::sqrt(oa * oa);
		const FT norm_ob = std::sqrt(ob * ob);
		assert(norm_oa != 0.0);
		assert(norm_ob != 0.0);
		return acos((oa * ob) / (norm_oa * norm_ob));
	}

    FT cotangent_from_h(Halfedge_handle hij) const
    {
        if (hij->is_border()) return 0.0;
        Vertex_handle vi = get_source_vertex(hij);
        Vertex_handle vj = get_target_vertex(hij);
        Vertex_handle vk = get_opposite_vertex(hij);
        return cotangent(vk->get_2d_point(),
                         vi->get_2d_point(),
                         vj->get_2d_point());
    }

	// at a
	FT cotangent(const Point_2& a, const Point_2& b, const Point_2& c) const
	{
		Vector_2 ab = b - a;
		Vector_2 ac = c - a;
		const FT cosine = ab*ac;
		const FT sine = determinant(ab, ac);
		return cosine/sine;
	}
    
	FT determinant(const Vector_2& a, const Vector_2& b) const
	{
		return a.x()*b.y() - a.y()*b.x();
	}
    
    FT cotangent_3d_from_h(Halfedge_handle hij) const
    {
        if (hij->is_border()) return 0.0;
        Vertex_handle vi = get_source_vertex(hij);
        Vertex_handle vj = get_target_vertex(hij);
        Vertex_handle vk = get_opposite_vertex(hij);
        return cotangent_3d(vk->get_3d_point(),
                            vi->get_3d_point(),
                            vj->get_3d_point());
    }
    
	// at a
	FT cotangent_3d(const Point_3& a, const Point_3& b, const Point_3& c) const
	{
		Vector_3 ab = b - a;
		Vector_3 ac = c - a;
        Vector_3 abxac = CGAL::cross_product(ab, ac);
		const FT cosine = ab*ac;
		const FT sine = std::sqrt(abxac * abxac);
		return cosine/sine;
	}
    
    // LOAD //
    
    FT compute_Fi(Vertex_handle vi) const
    {
        FT sum = 0.0;
        Halfedge_around_vertex_circulator hki = vi->vertex_begin();
		Halfedge_around_vertex_circulator end = hki;
		CGAL_For_all(hki, end)
        {
            if (hki->is_border()) continue;
            Facet_handle f = hki->facet();
            //sum += f->rho() * area_3d(f);
            sum += f->rho() * area(f);
        }
        return sum / 3.0;
    }
    
    Vector_3 compute_dFidxk(Halfedge_handle hki) const
    {
        Vector_3 sum_vec(0.0,0.0,0.0);
        
        if (!hki->is_border())
        {
            Facet_handle f = hki->facet();
            Vector_3 n = get_3d_normal(f);
            Vector_3 eij = get_3d_vector(hki->next());
            Vector_3 da = CGAL::cross_product(n, eij) / 6.0;
            sum_vec = sum_vec + f->rho() * da;
        }
        
        if (!hki->opposite()->is_border())
        {
            Facet_handle f = hki->opposite()->facet();
            Vector_3 n = get_3d_normal(f);
            Vector_3 eji = get_3d_vector(hki->opposite()->prev());
            Vector_3 da = CGAL::cross_product(n, eji) / 6.0;
            sum_vec = sum_vec + f->rho() * da;
        }
        
        return sum_vec;
    }
    
    Vector_3 compute_dFidxi(Vertex_handle vi) const
    {
        Vector_3 sum_vec(0.0,0.0,0.0);
        Halfedge_around_vertex_circulator hki = vi->vertex_begin();
        Halfedge_around_vertex_circulator end = hki;
        CGAL_For_all(hki, end)
        {
            if (hki->is_border()) continue;
            Facet_handle f = hki->facet();
            Vector_3 n = get_3d_normal(f);
            Vector_3 ejk = get_3d_vector(hki->prev());
            Vector_3 da = CGAL::cross_product(n, ejk) / 6.0;
            sum_vec = sum_vec + f->rho() * da;
        }
        return sum_vec;
    }
    
    // LIFTING //
    
	Point_3 lift_point(Facet_handle f, const Point_2& q) const
	{
		FT alpha, beta;
		get_barycentric_coordinates(q, f, alpha, beta);
		FT gamma = 1.0 - alpha - beta;
        
		Halfedge_handle he = f->halfedge();
		const Point_3 a = he->vertex()->get_3d_point();
		const Point_3 b = he->next()->vertex()->get_3d_point();
		const Point_3 c = he->next()->next()->vertex()->get_3d_point();
        
		return CGAL::ORIGIN + (alpha * (a - CGAL::ORIGIN) +
                               beta  * (b - CGAL::ORIGIN) +
                               gamma * (c - CGAL::ORIGIN) );
	}
    
	void get_barycentric_coordinates(const Point_2& query,
                                     Facet_handle f,
                                     FT& alpha,
                                     FT& beta) const
	{
		Halfedge_handle he = f->halfedge();
		const Point_2 p0 = he->vertex()->get_2d_point();
		const Point_2 p1 = he->next()->vertex()->get_2d_point();
		const Point_2 p2 = he->next()->next()->vertex()->get_2d_point();
        
		Triangle_2 triangle(p0, p1, p2);
		Triangle_2 q12(query, p1, p2);
		Triangle_2 q20(query, p2, p0);
        
		FT area = triangle.area();
		alpha = q12.area() / area;
		beta  = q20.area() / area;
	}

    // # ELEMENTS //
    
	unsigned nb_edges()
	{
		return this->size_of_halfedges()/2;
	}
    
    unsigned nb_vertices()
    {
        return this->size_of_vertices();
    }
    
    unsigned nb_facets()
    {
        return this->size_of_facets();
    }
    
    unsigned nb_grounded_vertices()
	{
		unsigned nb = 0;
		for (Vertex_iterator
             v = this->vertices_begin();
             v != this->vertices_end();
             v++)
		{
			if (is_vertex_grounded(v)) nb++;
		}
		return nb;
	}
    
    unsigned nb_free_vertices()
    {
        return nb_vertices() - nb_grounded_vertices();
    }
    
	unsigned nb_grounded_edges()
	{
		unsigned nb = 0;
		for (Edge_iterator
             he = this->edges_begin();
             he != this->edges_end();
             ++he)
		{
			if (is_edge_grounded(he)) nb++;
		}
		return nb;
	}
    
    unsigned nb_free_edges()
    {
        return nb_edges() - nb_grounded_edges();
    }
    
    unsigned nb_free_boundaries()
    {
		unsigned nb = 0;
		for (Edge_iterator
             he = this->edges_begin();
             he != this->edges_end();
             ++he)
		{
            if (is_free_boundary(he)) nb++;
		}
		return nb;
	}

    unsigned nb_borders()
    {
		unsigned nb = 0;
		for (Edge_iterator
             he = this->edges_begin();
             he != this->edges_end();
             ++he)
		{
            if (is_border_edge(he)) nb++;
		}
		return nb;
	}

    // TOPOLOGY COUNTING //

    void tag_vertices(const int tag)
	{
        for (Vertex_iterator
             v  = this->vertices_begin();
             v != this->vertices_end();
             v++)
        {
            v->tag() = tag;
		}
	}

    void set_all_vertex_boundary_indices(const int index)
	{
        for (Vertex_iterator
             v  = this->vertices_begin();
             v != this->vertices_end();
             v++)
        {
            v->boundary_index() = index;
		}
	}

	

    void constrain_vertices(const bool constrained)
	{
        for (Vertex_iterator
             v  = this->vertices_begin();
             v != this->vertices_end();
             v++)
        {
            v->constrained() = constrained;
		}
	}



	void tag_halfedges(const int tag)
	{
		for(Halfedge_iterator
			he = this->halfedges_begin();
			he != this->halfedges_end();
            he++)
		{
			he->tag() = tag;
		}
	}
    
	void tag_facets(const int tag)
	{
		for (Facet_iterator
             f = this->facets_begin();
             f != this->facets_end();
             f++)
		{
			f->tag() = tag;
		}
	}
    
	void tag_component(Facet_handle pSeedFacet,
                       const int tag_free,
                       const int tag_done)
	{
		pSeedFacet->tag() = tag_done;
		std::list<Facet_handle> facets;
		facets.push_front(pSeedFacet);
		while (!facets.empty())
		{
			Facet_handle pFacet = facets.front();
			facets.pop_front();
			pFacet->tag() = tag_done;
			Halfedge_around_facet_circulator pHalfedge = pFacet->facet_begin();
			Halfedge_around_facet_circulator end = pHalfedge;
			CGAL_For_all(pHalfedge,end)
			{
				Facet_handle pNFacet = pHalfedge->opposite()->facet();
				if (pNFacet != NULL && pNFacet->tag() == tag_free)
				{
					facets.push_front(pNFacet);
					pNFacet->tag() = tag_done;
				}
			}
		}
	}
    
	unsigned nb_components()
	{
		unsigned nb = 0;
		tag_facets(0);
		for (Facet_iterator
             pFacet = this->facets_begin();
             pFacet != this->facets_end();
             pFacet++)
		{
			if (pFacet->tag() == 0)
			{
				nb++;
				tag_component(pFacet,0,1);
			}
		}
		return nb;
	}
    
	unsigned nb_boundaries()
	{
		unsigned nb = 0;
		tag_halfedges(0);
		for (Halfedge_iterator
             he = this->halfedges_begin();
             he != this->halfedges_end();
             he++)
		{
			if (he->is_border() && he->tag() == 0)
			{
				nb++;
				Halfedge_handle curr = he;
				do
				{
					curr  = curr->next();
					curr->tag() = 1;
				}
				while(curr != he);
			}
		}
		return nb;
	}
    
	void trace_mesh()
	{
		std::cout << yellow << "INPUT PROPERTIES" << white << std::endl;
		std::cout << this->size_of_vertices() << " vertices" << std::endl;
		std::cout << this->size_of_facets() << " facets" << std::endl;
		std::cout << nb_edges() << " edges" << std::endl;
		std::cout << nb_boundaries() << " boundary(ies)" << std::endl;
		std::cout << nb_components() << " component(s)" << std::endl;
		trace_edge_len();
	}
    
	void trace_edge_len()
	{
		if (this->size_of_halfedges() == 0) return;
        
		int nb_edges = 0;
		FT sum_len = 0.0;
		FT max_len = 0.0;
		FT min_len = std::numeric_limits<double>::max();
		for (Edge_iterator he = this->edges_begin(); he != this->edges_end(); he++)
		{
			const FT edge_len = len_3d(he);
			sum_len += edge_len;
			min_len = std::min(min_len, edge_len);
			max_len = std::max(max_len, edge_len);
			nb_edges++;
		}
        
		std::cout << "Min edge length: " << min_len << std::endl;
		std::cout << "Max edge length: " << max_len << std::endl;
		std::cout << "Average edge length: " << sum_len / nb_edges << std::endl;
	}
    
    // CONVEX HULL //
    
	void update_2d_ch()
	{
		m_mesh_2d_ch.clear();
        
		std::vector<Point_2> points;
		for (Vertex_iterator
             v = this->vertices_begin();
             v != this->vertices_end(); v++)
		{
			points.push_back(v->get_2d_point());
		}
        
		CGAL::convex_hull_2(points.begin(),
                            points.end(),
                            std::back_inserter(m_mesh_2d_ch));
	}
    
	template <class OutputIterator>
	void get_2d_ch_facets_as_3d_triangles(OutputIterator out)
	{
		int size_ch = int(m_mesh_2d_ch.size());
		if (size_ch < 3) return;
        
		// select pivot vertex
		typename std::list<Point_2>::iterator it = m_mesh_2d_ch.begin();
		const Point_3& pivot = to_3d_no_height(*it);
		it++; // will start at second point
		for(int i=0; i<(size_ch-2); i++)
		{
			const Point_2& p = *it;
			it++; // move to next and stay there
			const Point_2& q = *it;
			const Triangle_3 triangle(pivot, to_3d_no_height(p), to_3d_no_height(q));
			*out++ = triangle;
		}
	}
    
    Point_3 to_3d_no_height(const Point_2& p) const
	{
		return Point_3(p.x(), p.y(), 0.0);
	}
    
	// DUAL //
    
	Point_2 get_face_cw(Facet_handle f) const
	{
		Halfedge_handle hij = f->halfedge();
		const Point_2 pi = get_source_vertex(hij)->get_2d_point();
		const Vector_2 eij = get_vector_2d(hij);
		const Vector_2 eij90 = rotate90(eij);
		const FT dij = get_dij(hij);
		const FT hk = get_hk(hij);
		const FT lij = len_2d(hij);
		return pi + (dij/lij)*eij + (hk/lij)*eij90;
	}
    
	Point_2 get_edge_cw(Halfedge_handle hij) const
	{
		Vertex_handle vi = get_source_vertex(hij);
		Vertex_handle vj = get_target_vertex(hij);

        Point_2 qi = vi->get_2d_point();
		Point_2 qj = vj->get_2d_point();

		FT lij = len_2d(hij);
		FT dij = get_dij(hij);
		return qi + (dij/lij)*(qj - qi);
	}
    
    FT get_omega(Halfedge_handle hij) const
    {
		Vertex_handle vi = get_source_vertex(hij);
		Vertex_handle vj = get_target_vertex(hij);
		FT wi = vi->weight();
		FT wj = vj->weight();
        FT omega = wj - wi;
        
        const std::vector<Form>& forms = hij->forms();
        for (unsigned i = 0; i < forms.size(); ++i)
        {
            // normal is integrable
            omega += m_alpha[i] * forms[i].t();
        }
        
        return omega;
    }
    
    Vector_2 get_face_vector(Facet_handle f) const
    {
        Halfedge_handle hij = f->halfedge();
        
        Vector_2 ei = get_vector_2d(hij->next());
        Vector_2 ej = get_vector_2d(hij->prev());
        
        Vector_2 rei = rotate90(ei);
        Vector_2 rej = rotate90(ej);
        
        FT omega_ki = get_omega(hij->prev());
        FT omega_jk = get_omega(hij->next());

        FT a = area(f);
        return (omega_ki*rei - omega_jk*rej) / (2.0*a);
    }
    
	FT get_dij(Halfedge_handle hij) const
	{
        FT lij = len_2d(hij);
		FT Wij = get_omega(hij);
		return 0.5*(lij - Wij/lij);
	}
    
	FT get_hk(Halfedge_handle hij) const
	{
		if (hij->is_border()) return 0.0;

		Vertex_handle vi = get_source_vertex(hij);
		Vertex_handle vj = get_target_vertex(hij);
		Vertex_handle vk = get_opposite_vertex(hij);

		Point_2 pi = vi->get_2d_point();
		Point_2 pj = vj->get_2d_point();
		Point_2 pk = vk->get_2d_point();
        
        FT Wkj = - get_omega(hij->next());
        FT Wki =   get_omega(hij->prev());
        
		const FT lij = len_2d(hij);
		const FT coti = cotangent(pi, pj, pk);
		const FT cotj = cotangent(pj, pk, pi);
		const FT cotk = cotangent(pk, pi, pj);
        
		return 0.5*(cotk*lij + Wki*cotj/lij + Wkj*coti/lij);
	}
    
    FT get_unweighted_hk(Halfedge_handle hij) const
    {
        if (hij->is_border()) return 0.0;

		Vertex_handle vi = get_source_vertex(hij);
		Vertex_handle vj = get_target_vertex(hij);
		Vertex_handle vk = get_opposite_vertex(hij);

		const Point_2 pi = vi->get_2d_point();
		const Point_2 pj = vj->get_2d_point();
		const Point_2 pk = vk->get_2d_point();
        
		const FT lij = len_2d(hij);
		const FT cotk = cotangent(pk, pi, pj);
		return 0.5*cotk*lij;
	}
    
    FT get_surf_hk(Halfedge_handle hij) const
    {
        if (hij->is_border()) return 0.0;
        
		Vertex_handle vi = get_source_vertex(hij);
		Vertex_handle vj = get_target_vertex(hij);
		Vertex_handle vk = get_opposite_vertex(hij);
        
		const Point_3 pi = vi->get_3d_point();
		const Point_3 pj = vj->get_3d_point();
		const Point_3 pk = vk->get_3d_point();
        
		const FT lij = len_3d(hij);
		const FT cotk = cotangent_3d(pk, pi, pj);
		return 0.5*cotk*lij;
	}
    
    FT get_star1(Halfedge_handle hij) const
	{
		FT dual_len = get_hk(hij) + get_hk(hij->opposite());
		FT primal_len = len_2d(hij);
		FT star = dual_len / primal_len;
        return star;
	}
    
    FT get_unweighted_star1(Halfedge_handle hij) const
	{
		FT dual_len = get_unweighted_hk(hij) + get_unweighted_hk(hij->opposite());
		FT primal_len = len_2d(hij);
		FT star = dual_len / primal_len;
        return star;
	}

    FT get_surf_star1(Halfedge_handle hij) const
	{
		FT dual_len = get_surf_hk(hij) + get_surf_hk(hij->opposite());
		FT primal_len = len_3d(hij);
		FT star = dual_len / primal_len;
        return star;
	}

	FT get_star0(Vertex_handle vi) const
	{
		FT sum = 0.0;
		const FT wi = vi->weight();
		Halfedge_around_vertex_circulator hji = vi->vertex_begin();
		Halfedge_around_vertex_circulator end = hji;
		CGAL_For_all(hji, end)
		{
            FT dual_len = get_hk(hji) + get_hk(hji->opposite());
            FT dij = get_dij(hji->opposite());
            sum += 0.5*dij*dual_len;
        }
        return sum;
	}

    FT get_unweighted_star0(Vertex_handle vi) const
	{
		FT sum = 0.0;
		const FT wi = vi->weight();
		Halfedge_around_vertex_circulator hji = vi->vertex_begin();
		Halfedge_around_vertex_circulator end = hji;
		CGAL_For_all(hji, end)
		{
            FT dual_len = get_unweighted_hk(hji) + get_unweighted_hk(hji->opposite());
            FT dij = 0.5*len_2d(hji);
            sum += 0.5*dij*dual_len;
        }
        return sum;
	}

	// primal x dual > 0
	Segment_2 get_unbounded_dual_edge(Halfedge_handle hij, bool with_bdry) const
	{
		Halfedge_handle hji = hij->opposite();
		bool left_inside  = !hij->is_border();
		bool right_inside = !hji->is_border();
        
        if (with_bdry && left_inside != right_inside)
        {
            Point_2 left_cw, right_cw;
            if (left_inside) left_cw = get_face_cw(hij->facet());
            if (right_inside) right_cw = get_face_cw(hji->facet());
            
            FT step = get_star1(hij) + hij->sigma();
            Vector_2 eij = get_vector_2d(hij);
            Vector_2 eij90 = rotate90(eij);
            
            if (!left_inside) left_cw = right_cw + step*eij90;
            if (!right_inside) right_cw = left_cw - step*eij90;
            return Segment_2(right_cw, left_cw);
        }
        
		if (!left_inside || !right_inside)
			return Segment_2(CGAL::ORIGIN, CGAL::ORIGIN);
        
        Point_2 left_cw  = get_face_cw(hij->facet());
        Point_2 right_cw = get_face_cw(hji->facet());
        return Segment_2(right_cw, left_cw);
	}
    
	// primal x dual > 0
	Segment_2 get_bounded_dual_edge(Halfedge_handle hij) const
	{
		Point_2 left_cw;
		if (hij->is_border()) left_cw = get_edge_cw(hij);
		else left_cw = get_face_cw(hij->facet());
        
		Point_2 right_cw;
		Halfedge_handle hji = hij->opposite();
		if (hji->is_border()) right_cw = get_edge_cw(hji);
		else right_cw = get_face_cw(hji->facet());
        
		return Segment_2(right_cw, left_cw);
	}
    
	void get_dual_cell(Vertex_handle v, std::vector<Point_2>& cell) const
	{
		// traverse in clockwise
		Halfedge_around_vertex_circulator he = v->vertex_begin();
		Halfedge_around_vertex_circulator end = he;
		CGAL_For_all(he, end)
		{
			if (he->opposite()->is_border())
			{
				cell.push_back(v->get_2d_point());
				cell.push_back(get_edge_cw(he));
			}
            
			if (he->is_border())
				cell.push_back(get_edge_cw(he));
			else
				cell.push_back(get_face_cw(he->facet()));
		}
		std::reverse(cell.begin(), cell.end());
	}
    
    void get_barycentric_cell(Vertex_handle v, std::vector<Point_2>& cell) const
	{
		// traverse in clockwise
		Halfedge_around_vertex_circulator he = v->vertex_begin();
		Halfedge_around_vertex_circulator end = he;
		CGAL_For_all(he, end)
		{
			if (he->is_border()) continue;
            cell.push_back(CGAL::centroid(triangle(he->facet())));
		}
		std::reverse(cell.begin(), cell.end());
	}
    
	// TOPOLOGICAL OPERATORS //
    
	bool has_edge(Halfedge_handle heq)
	{
		for (Halfedge_iterator
             he = this->halfedges_begin();
             he != this->halfedges_end();
             he++)
		{
			if (he == heq) return true;
		}
		return false;
	}
    
	bool has_vertex(Vertex_handle vq)
	{
		for (Vertex_iterator
             v = this->vertices_begin();
             v != this->vertices_end();
             v++)
		{
			if (v == vq) return true;
		}
		return false;
	}
    
	bool join_facets_before_collapse(Halfedge_handle he)
	{
		// precondition: he is not on border
		assert(!he->is_border());
        
		if (he->opposite()->is_border())
			return join_facets_before_collapse_boundary_case(he);
		else
			return join_facets_before_collapse_inner_case(he);
	}
    
	bool valid_join_facet(Halfedge_handle he)
	{
		if(circulator_size(he->vertex_begin()) < 3) return false;
		if(circulator_size(he->opposite()->vertex_begin()) < 3) return false;
		return true;
	}
    
	bool valid_join_facet_pair(Facet_handle f1, Facet_handle f2)
	{
		if(f1 == f2) return false;
		if(f1 == Facet_handle()) return false;
		if(f2 == Facet_handle()) return false;
		return true;
	}
    
	bool join_facets_before_collapse_inner_case(Halfedge_handle he)
	{
		Halfedge_handle qr = he->next();
		Halfedge_handle rp = he->prev();
		Halfedge_handle ps = he->opposite()->next();
		Halfedge_handle sq = he->opposite()->prev();
        
		Facet_handle frq = qr->opposite()->facet();
		Facet_handle fpr = rp->opposite()->facet();
		Facet_handle fsp = ps->opposite()->facet();
		Facet_handle fqs = sq->opposite()->facet();
        
		// check if faces are different
		if(valid_join_facet_pair(fsp, frq))
		{
			if(valid_join_facet(ps) && valid_join_facet(qr))
			{
				this->join_facet(ps);
				this->join_facet(qr);
				return true;
			}
		}
        
		if(valid_join_facet_pair(fsp, fpr))
		{
			if(valid_join_facet(ps) && valid_join_facet(rp))
			{
				this->join_facet(ps);
				this->join_facet(rp);
				return true;
			}
		}
        
		if(valid_join_facet_pair(fqs, frq))
		{
			if(valid_join_facet(sq) && valid_join_facet(qr))
			{
				this->join_facet(sq);
				this->join_facet(qr);
				return true;
			}
		}
        
		if(valid_join_facet_pair(fqs, fpr))
		{
			if(valid_join_facet(sq) && valid_join_facet(rp))
			{
				this->join_facet(sq);
				this->join_facet(rp);
				return true;
			}
		}
        
		std::cerr << red << "join facets failed" << white << std::endl;
		return false;
	}
    
	// return true if relocating v to new_pos would
	// preserve embedding
	bool valid_embedding(Vertex_handle v, const Point_2& new_pos)
	{
		Halfedge_around_vertex_circulator he = v->vertex_begin();
		Halfedge_around_vertex_circulator end = he;
		CGAL_For_all(he, end)
		{
			// skip boundary edges
			if (he->is_border())
				continue;
            
			// get three points
			const Point_2 p = he->vertex()->get_2d_point();
			const Point_2 q = he->next()->vertex()->get_2d_point();
			const Point_2 r = he->next()->next()->vertex()->get_2d_point();
			Triangle_2 triangle(p, q, r);
			CGAL::Orientation orientation = triangle.orientation();
            
			Triangle_2 new_triangle(new_pos, q, r);
			CGAL::Orientation new_orientation = new_triangle.orientation();
            
			if(new_orientation != orientation)
				return false;
		}
		return true;
	} // end valid_embedding
    
    
	bool join_facets_before_collapse_boundary_case(Halfedge_handle he)
	{
		// given edge = pq
		CGAL_precondition(he->opposite()->is_border());
        
		Halfedge_handle qr = he->next();
		Halfedge_handle rp = he->next()->next();
        
		Facet_handle frq = qr->opposite()->facet();
		Facet_handle fpr = rp->opposite()->facet();
        
		if(frq != Facet_handle())
		{
			if(valid_join_facet(qr))
			{
				this->join_facet(qr);
				return true;
			}
		}
		else if(fpr != Facet_handle())
		{
			if(valid_join_facet(rp))
			{
				this->join_facet(rp);
				return true;
			}
		}
		std::cerr << red << "never come here: in join_facets_before_collapse_boundary_case" << std::endl;
		return false;
	}
    
	bool is_regular_flippable(Halfedge_handle he) // edge pq
	{
		// do not flip boundary edges
		if (is_border_edge(he)) return false;
        
		// do not flip edge if the flipped edge already exist
		Vertex_handle p = he->vertex();
		Vertex_handle q = he->opposite()->vertex();
		Vertex_handle r = he->next()->vertex();
		Vertex_handle s = he->opposite()->next()->vertex();
		if (are_neighbors(r,s)) return false;
        
		return edge_flip_would_improve_regular(he);
	}
    
	// edge ac flippable if ...
	bool edge_flip_would_improve_regular(Halfedge_handle he)
	{
		Vertex_handle va = he->vertex();
		Vertex_handle vb = he->next()->vertex();
		Vertex_handle vc = he->opposite()->vertex();
		Vertex_handle vd = he->opposite()->next()->vertex();
        
		Point_2 a = va->get_2d_point();
		Point_2 b = vb->get_2d_point();
		Point_2 c = vc->get_2d_point();
		Point_2 d = vd->get_2d_point();
        
		// before flip
		Triangle_2 abc(a, b, c);
		Triangle_2 acd(a, c, d);
        
		// after flip
		Triangle_2 bcd(b, c, d);
		Triangle_2 abd(a, b, d);
        
		// check preservation of orientation (pairwise)
		if(abc.orientation() != bcd.orientation()) return false;
		if(abc.orientation() != abd.orientation()) return false;
		if(acd.orientation() != bcd.orientation()) return false;
		if(acd.orientation() != abd.orientation()) return false;
        
		// check power of b with respect to face acd
		// see line 1978 of Regular_triangulation_2.h of CGAL 4.0
		if (power_test(va, vc, vd,  vb) != CGAL::ON_POSITIVE_SIDE)
			return false;
        return true;
	}
    
	// check power of vertex vq with respect to face (v0-v1-v2)
	CGAL::Oriented_side power_test(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Vertex_handle vq)
	{
		Point_2 p0 = v0->get_2d_point(); FT w0 = v0->weight();
		Point_2 p1 = v1->get_2d_point(); FT w1 = v1->weight();
		Point_2 p2 = v2->get_2d_point(); FT w2 = v2->weight();
		Point_2 pq = vq->get_2d_point(); FT wq = vq->weight();
		return CGAL::power_testC2(p0.x(), p0.y(), w0,
                                  p1.x(), p1.y(), w1,
                                  p2.x(), p2.y(), w2,
                                  pq.x(), pq.y(), wq);
	}
    
	// return true if vertex v is submerged
	bool is_submerged(Vertex_handle v)
	{
		Halfedge_around_vertex_circulator he = v->vertex_begin();
		Halfedge_around_vertex_circulator end = he;
		CGAL_For_all(he, end)
        if(he->wstar() >= 0.0)
            return false; // all edges must have negative
		// stars to qualify the vertex as submerged
		return true;
	}

	bool is_intrinsic_delaunay_flippable(Halfedge_handle he) // edge pq
	{
		// do not flip boundary edges
		if (is_border_edge(he)) return false;
        
		// do not flip edge if the flipped edge already exist
		Vertex_handle p = he->vertex();
		Vertex_handle q = he->opposite()->vertex();
		Vertex_handle r = he->next()->vertex();
		Vertex_handle s = he->opposite()->next()->vertex();
		if (are_neighbors(r,s)) return false;
        
		return intrinsic_edge_flip_would_improve_angles(he);
	}
    
	// edge ac flippable if increases the minimum angle of triangles
	// adjacent to the flipped edge.
	// check orientations too
	bool intrinsic_edge_flip_would_improve_angles(Halfedge_handle he)
	{
		const Point_3& a = he->vertex()->point();
		const Point_3& b = he->next()->vertex()->point();
		const Point_3& c = he->opposite()->vertex()->point();
		const Point_3& d = he->opposite()->next()->vertex()->point();
        
		// check negative area
		Triangle_3 abc(a, b, c);
		Triangle_3 acd(a, c, d);
        
		Triangle_3 bcd(b, c, d);
		Triangle_3 abd(a, b, d);
        
		// check preservation of orientation (pairwise)
		//if(abc.orientation() != bcd.orientation()) return false;
		//if(abc.orientation() != abd.orientation()) return false;
		//if(acd.orientation() != bcd.orientation()) return false;
		//if(acd.orientation() != abd.orientation()) return false;

		FT min_angle_before = compute_min_angle_3d(a, b, c, d);
		FT min_angle_after  = compute_min_angle_3d(d, a, b, c);

		FT diff = min_angle_after - min_angle_before;
        
		return diff > 1e-3; // FIXME should be zero but there are constructions...
	}

	// compute increase of minimum angle of triangles incident to he, upon flip.
	// if the angle decreases then the function returns 0.
	double angle_increase_upon_flip(Halfedge_handle he)
	{
		const Point_3& a = he->vertex()->point();
		const Point_3& b = he->next()->vertex()->point();
		const Point_3& c = he->opposite()->vertex()->point();
		const Point_3& d = he->opposite()->next()->vertex()->point();
        
		// check negative area
		Triangle_3 abc(a, b, c);
		Triangle_3 acd(a, c, d);
        
		Triangle_3 bcd(b, c, d);
		Triangle_3 abd(a, b, d);
        
		FT min_angle_before = compute_min_angle_3d(a, b, c, d);
		FT min_angle_after  = compute_min_angle_3d(d, a, b, c);

		FT increase = min_angle_after - min_angle_before;
		return increase < 0.0 ? 0.0 : increase;
	}


	bool is_delaunay_flippable(Halfedge_handle he) // edge pq
	{
		// do not flip boundary edges
		if (is_border_edge(he)) return false;
        
		// do not flip edge if the flipped edge already exist
		Vertex_handle p = he->vertex();
		Vertex_handle q = he->opposite()->vertex();
		Vertex_handle r = he->next()->vertex();
		Vertex_handle s = he->opposite()->next()->vertex();
		if (are_neighbors(r,s)) return false;
        
		return edge_flip_would_improve_angles(he);
	}
    
	// edge ac flippable if increases the minimum angle of triangles
	// adjacent to the flipped edge.
	// check negative areas too
	bool edge_flip_would_improve_angles(Halfedge_handle he)
	{
		Point_2 a = he->vertex()->get_2d_point();
		Point_2 b = he->next()->vertex()->get_2d_point();
		Point_2 c = he->opposite()->vertex()->get_2d_point();
		Point_2 d = he->opposite()->next()->vertex()->get_2d_point();
        
		// check negative area
		Triangle_2 abc(a, b, c);
		Triangle_2 acd(a, c, d);
        
		Triangle_2 bcd(b, c, d);
		Triangle_2 abd(a, b, d);
        
		// check preservation of orientation (pairwise)
		if(abc.orientation() != bcd.orientation()) return false;
		if(abc.orientation() != abd.orientation()) return false;
		if(acd.orientation() != bcd.orientation()) return false;
		if(acd.orientation() != abd.orientation()) return false;
        
		if(CGAL::side_of_bounded_circle(a,b,c, d) == CGAL::ON_BOUNDED_SIDE)
			return true;
		if(CGAL::side_of_bounded_circle(d,a,c, b) == CGAL::ON_BOUNDED_SIDE)
			return true;
        
		return false;
	}
    
	// returns true if boundary edge can be split at new_pos
	// while preserving a valid embedding.
	bool valid_embedding_upon_boundary_split(Halfedge_handle he, const Point_2& new_pos)
	{
		// pick inner half edge
		CGAL_assertion(he->is_border() || he->opposite()->is_border());
        
		Halfedge_handle h = he->is_border() ? he->opposite() : he;
        
		Vertex_handle vp = get_source_vertex(h);
		Vertex_handle vq = get_target_vertex(h);
		Vertex_handle vr = get_opposite_vertex(h);
        
		const Point_2 p = vp->get_2d_point();
		const Point_2 q = vq->get_2d_point();
		const Point_2 r = vr->get_2d_point();
        
		// valid  = same orientation for two 2D triangles upon edge split
		Triangle_2 nqr(new_pos, q, r);
		Triangle_2 nrp(new_pos, r, p);
		return nqr.orientation() == nrp.orientation();
	}
    
	// check if collapsing edge he will preserve valid embedding
	bool valid_embedding_upon_collapse(Halfedge_handle he, const Point_2& new_pos)
	{
		Vertex_handle vp = get_source_vertex(he);
		Vertex_handle vq = get_target_vertex(he);
		if(!valid_incident_facet_embedding_upon_relocation(vp, vq, new_pos))
			return false;
		if(!valid_incident_facet_embedding_upon_relocation(vq, vp, new_pos))
			return false;
		return true;
	}
    
	// check if all faces incident to vp (but the one containing vx) preserve
	// orientation when vp is changed to new pos.
	// note that the (one or two) faces incident to both vp and vx
	// gets deleted upon collapse.
	bool valid_incident_facet_embedding_upon_relocation(Vertex_handle vp,
                                                        Vertex_handle vx, const Point_2& new_pos)
	{
		Halfedge_around_vertex_circulator he = vp->vertex_begin();
		Halfedge_around_vertex_circulator end = he;
		CGAL_For_all(he, end)
		{
			if(he->is_border())
				continue;
            
			const Point_2 p = he->vertex()->get_2d_point();
            
			// skip face containing vx
			Vertex_handle vq = he->next()->vertex();
			Vertex_handle vr = he->next()->next()->vertex();
			if(vq == vx || vr == vx)
				continue;
            
			const Point_2 q = vq->get_2d_point();
			const Point_2 r = vr->get_2d_point();
			Triangle_2 pqr(p, q, r);
			Triangle_2 new_pqr(new_pos, q, r);
			if(pqr.orientation() != new_pqr.orientation())
				return false;
		}
		return true;
	}

	// assume edge ac
	// compute min angle of incident triangles
	// b,d are two points facing ac
	FT compute_min_angle_3d(const Point_3& a, const Point_3& b,
                            const Point_3& c, const Point_3& d)
	{
		FT cab = angle_3d(c, a , b);
		FT bca = angle_3d(b, c , a);
		FT abc = angle_3d(a, b , c);
        
		FT dac = angle_3d(d, a , c);
		FT cda = angle_3d(c, d , a);
		FT acd = angle_3d(a, c , d);
        
		return std::min(std::min(cab, std::min(bca, abc)),
                        std::min(dac, std::min(cda, acd)));
	}
    
	// assume edge ac
	// compute min angle of incident triangles
	// b,d are two points facing ac
	FT compute_min_angle_2d(const Point_2& a, const Point_2& b,
                            const Point_2& c, const Point_2& d)
	{
		FT cab = angle_2d(c, a , b);
		FT bca = angle_2d(b, c , a);
		FT abc = angle_2d(a, b , c);
        
		FT dac = angle_2d(d, a , c);
		FT cda = angle_2d(c, d , a);
		FT acd = angle_2d(a, c , d);
        
		return std::min(std::min(cab, std::min(bca, abc)),
                        std::min(dac, std::min(cda, acd)));
	}
    
	bool is_collapsible(Halfedge_handle he) const
	{
		// precondition: h is not a boundary halfedge
		// (its opposite halfedge can be)
        
		if (he->is_border())
		{
			std::cerr << red << "boundary edge given to collapse." << white << std::endl;
			return false;
		}
        
		Halfedge_handle ho = he->opposite();
        
		// opposite halfedge is on boundary
		if (ho->is_border())
		{
			Vertex_handle vr = get_opposite_vertex(he);
			Vertex_handle vs = get_target_vertex(ho->next());
			Vertex_handle vt = get_source_vertex(ho->prev());
            
			if (vs == vt) return false; // do not close degree-3 hole
			if (!check_link_test(he)) return false;
			return true;
		}
        
		// inner case
		Vertex_handle vr = get_opposite_vertex(he);
		Vertex_handle vs = get_opposite_vertex(ho);
		if (vr == vs) return false;
        
		// inner edge but two boundary vertices
		if (is_border_vertex(he->vertex()) &&
			is_border_vertex(ho->vertex()))
			return false;
        
		if (!check_link_test(he)) return false;
        
		return true;
	}
    
	bool is_border_edge(Halfedge_handle he) const
	{
		if (he->is_border() || he->opposite()->is_border()) return true;
		return false;
	}
    
	bool is_border_vertex(Vertex_handle v) const
	{
		Halfedge_around_vertex_circulator vcirc = v->vertex_begin();
		Halfedge_around_vertex_circulator vend = vcirc;
		CGAL_For_all(vcirc, vend)
		{
			if (vcirc->is_border()) return true;
		}
		return false;
	}
    
	// return true if edge has at least one vertex on border
	bool touch_border(Halfedge_handle he)
	{
		if(is_border_edge(he)) return true;
		if(is_border_vertex(get_source_vertex(he))) return true;
		if(is_border_vertex(get_target_vertex(he))) return true;
		return false;
	}
    
	bool check_link_test(Halfedge_handle h) const
	{
		// h = (p -> q)
		Vertex_handle vp = get_source_vertex(h);
		Vertex_handle vq = get_target_vertex(h);
        
		if(h->opposite()->is_border())
		{
			// (p, q, r)
			Vertex_handle vr = get_opposite_vertex(h);
			Halfedge_around_vertex_circulator vcirc = vp->vertex_begin();
			Halfedge_around_vertex_circulator vend = vcirc;
			CGAL_For_all(vcirc, vend)
			{
				Vertex_handle v = get_source_vertex(vcirc);
				if(v == vq || v == vr) continue;
				if(are_neighbors(v, vq)) return false;
			}
			return true;
		}
        
		// (p, q, r)
		Vertex_handle vr = get_opposite_vertex(h);
		// (q, p, s)
		Vertex_handle vs = get_opposite_vertex(h->opposite());
		//
		Halfedge_around_vertex_circulator vcirc = vp->vertex_begin();
		Halfedge_around_vertex_circulator vend  = vcirc;
		CGAL_For_all(vcirc, vend)
		{
			Vertex_handle v = get_source_vertex(vcirc);
			if (v == vq || v == vr || v == vs) continue;
			if (are_neighbors(v, vq)) return false;
		}
		return true;
	}
    
	bool are_neighbors(Vertex_handle va, Vertex_handle vb) const
	{
		Halfedge_around_vertex_circulator vcirc = va->vertex_begin();
		Halfedge_around_vertex_circulator vend  = vcirc;
		CGAL_For_all(vcirc, vend)
		{
			if (get_source_vertex(vcirc) == vb)
				return true;
		}
		return false;
	}
    
	// IO //
    
	template <class OutputIterator>
	void get_3d_points(OutputIterator out)
	{
		for (Vertex_iterator v = this->vertices_begin(); v != this->vertices_end(); v++)
			*out++ = v->point();
	}
    
	template <class OutputIterator>
	void get_boundary_2d_points(OutputIterator out)
	{
		for (Vertex_iterator v = this->vertices_begin(); v != this->vertices_end(); v++)
			if(is_border_vertex(v))
				*out++ = v->get_2d_point(); // add 2D point
	}
    
	// get all 2D points
	template <class OutputIterator>
	void get_2d_points(OutputIterator out)
	{
		for (Vertex_iterator v = this->vertices_begin(); v != this->vertices_end(); v++)
			*out++ = v->get_2d_point(); // add 2D point
	}
    
	template <class EdgeOutputIterator, class FacetOutputIterator>
	void get_edge_and_facet_indices(EdgeOutputIterator edge_out,
                                    FacetOutputIterator facet_out)
	{
		// set vertex indices
		set_vertex_indices();
        
		// edges
		for (Edge_iterator he = this->edges_begin(); he != this->edges_end(); he++)
		{
			*edge_out++ = get_source_vertex(he)->index();
			*edge_out++ = get_target_vertex(he)->index();
		}
        
		// facets
		for (Facet_iterator f = this->facets_begin(); f != this->facets_end(); f++)
		{
			Halfedge_around_facet_circulator he = f->facet_begin();
			Halfedge_around_facet_circulator end = he;
			CGAL_For_all(he, end)
            *facet_out++ = get_target_vertex(he)->index();
		}
	}

	template <class OutputIterator>
	void get_projected_facets(OutputIterator out, const float opacity)
	{
		typedef typename Render_transparent<Kernel>::Projected_facet Proj_face;
		for (Facet_iterator f = this->facets_begin(); f != this->facets_end(); f++)
			*out++ = Proj_face(triangle_3d(f), opacity);
	}
    
	template <class OutputIterator>
	void get_boundary_edge_indices(OutputIterator out)
	{
		set_vertex_indices();
		for (Edge_iterator he = this->edges_begin(); he != this->edges_end(); he++)
		{
			// careful: edge iterator provides only 1 halfedge out of 2 (arbitrarily)
			if(he->is_border() || he->opposite()->is_border())
			{
				*out++ = get_source_vertex(he)->index();
				*out++ = get_target_vertex(he)->index();
			}
		}
	}
    
	void set_vertex_indices()
	{
		int index = 0;
		for (Vertex_iterator v = this->vertices_begin(); v != this->vertices_end(); v++)
			v->index() = index++;
	}
    
	// Height probing //
    
	Triangle_2 get_half_2d_diamond(Halfedge_handle he) const
	{
		CGAL_assertion(!he->is_border());
		const Point_2 p = get_source_vertex(he)->get_2d_point();
		const Point_2 q = get_target_vertex(he)->get_2d_point();
		const Point_2 b = CGAL::centroid(triangle(he->facet()));
		return Triangle_2(p, q, b);
	}
    
	FT evaluate_height(Facet_handle f, const Point_2& query) const
	{
		FT beta = 0.0;
		FT alpha = 0.0;
		get_barycentric_coordinates(query, f, alpha, beta);
		FT gamma = 1.0 - alpha - beta;
        
		Halfedge_handle he = f->halfedge();
		Vertex_handle v0 = he->vertex();
		Vertex_handle v1 = he->next()->vertex();
		Vertex_handle v2 = he->next()->next()->vertex();
		return (alpha*v0->height() + beta*v1->height() + gamma*v2->height());
	}
    
	template <class Oracle>
    FT evaluate_height_error(Halfedge_handle he, Oracle* oracle, const unsigned nb_probes) const
    {
        FT E0 = evaluate_height_error_per_halfedge(he, oracle, nb_probes);
        FT E1 = evaluate_height_error_per_halfedge(he->opposite(), oracle, nb_probes);
        return E0 + E1;
    }
    
    template <class Oracle>
    FT evaluate_height_error_per_halfedge(Halfedge_handle he, Oracle* oracle, const unsigned nb_probes) const
    {
        if (he->is_border()) return 0.0;
        
		typedef std::vector<Point_2> Point_list;
		typedef std::back_insert_iterator<Point_list> OutputIterator;
        
		std::vector<Point_2> samples;
		Triangle_2 triangle = get_half_2d_diamond(he);
		random_uniform_sample<Kernel,OutputIterator>(triangle,
                                                     nb_probes,
                                                     std::back_inserter(samples));
        
		FT sum_diff = 0.0;
		unsigned inside_samples = 0;
		Facet_handle f = he->facet();
		for (unsigned i = 0; i < samples.size(); ++i)
		{
			bool inside = false;
			const Point_2& q = samples[i];
			const FT h1 = oracle->evaluate_height_check(q, inside);
			if (!inside) continue;
            
            const FT h2 = evaluate_height(f, q);
            const FT diff = std::abs(h1 - h2);
            sum_diff += diff;
            inside_samples++;
		}
        
		if (inside_samples == 0) return 0.0;
        return sum_diff / FT(inside_samples);
    }
    
	// assuming he is the inner halfedge of a boundary edge
	// return true if angle facing he is not acute
	bool boundary_edge_to_conform(Halfedge_handle he)
	{
		// skip border edge
		if(!is_border_edge(he)) return false;
        
		// pick inner halfedge of boundary edge
		Halfedge_handle h = he->is_border() ? he->opposite() : he;
		return !is_2d_acute_angle(h->next());
	}
    
	bool is_2d_acute_angle(Halfedge_handle he)
	{
		const Point_2& p = he->vertex()->get_2d_point();
		const Point_2& q = he->prev()->vertex()->get_2d_point();
		const Point_2& r = he->next()->vertex()->get_2d_point();
		return (q-p) * (r-p) > 0.0;
	}
    
	template <class Oracle>
	bool must_be_refined_for_height(Halfedge_handle he, Oracle* oracle,
                                    const unsigned nb_probes, const double max_height) const
	{
		FT E = evaluate_height_error(he, oracle, nb_probes, max_height);
        return (E > max_height);
	}
    
	template <class Oracle>
	void debug_draw_edges_probing(Oracle* oracle,
                                  const unsigned nb_probes,
                                  const FT height_tolerance)
	{
		srand(0);
		glDisable(GL_LIGHTING);
		glColor3ub(255,0,0);
        
		glBegin(GL_LINES);
		for (Halfedge_iterator
             he  = this->halfedges_begin();
             he != this->halfedges_end();
             he++)
		{
			if (he->is_border()) continue;
			if (must_be_refined_for_height(he, oracle, nb_probes, height_tolerance))
			{
				const Point_3& p = he->vertex()->point();
				const Point_3& q = he->opposite()->vertex()->point();
				glVertex3d(p.x(), p.y(), p.z());
				glVertex3d(q.x(), q.y(), q.z());
			}
		}
		glEnd();
	}
    
	template <class Oracle>
	void draw_edge_probing(Oracle* oracle,
                           const unsigned nb_probes,
                           Halfedge_handle he)
	{
		typedef std::vector<Point_2> Point_list;
		typedef std::back_insert_iterator<Point_list> OutputIterator;
        
		// probe half-diamond triangle
		Triangle_2 triangle = get_half_2d_diamond(he);
		std::vector<Point_2> samples;
		random_uniform_sample<Kernel,OutputIterator>(triangle,
                                                     nb_probes,
                                                     std::back_inserter(samples));
        
		// draw line segments to depict the height diff.
		glDisable(GL_LIGHTING);
		glBegin(GL_LINES);
		Facet_handle f = he->facet();
		for (unsigned i = 0; i < samples.size(); ++i)
		{
			const Point_2& query = samples[i];
			const FT h1 = oracle->evaluate_height(query);
			const FT h2 = evaluate_height(f, query);
            
			const float diff = std::fabs(h1 - h2);
			glColor3f(1.0f - diff, 0.0, diff);
            
			glVertex3d(query.x(), query.y(), h1);
			glVertex3d(query.x(), query.y(), h2);
		}
		glEnd();
	}
	
    // HARMONIC FORMS //
	
	void reset_forms()
	{
        m_alpha.clear();
        
		for (Halfedge_iterator he = this->halfedges_begin(); he != this->halfedges_end(); he++)
            he->forms().clear();

		for (Vertex_iterator v = this->vertices_begin(); v != this->vertices_end(); v++)
			v->forms().clear();
	}
    
    void init_alpha(int nb)
    {
        m_alpha = std::vector<FT>(nb, 0.0);
    }

	// MOUSE PICKING

	bool toggle_boundary_vertex_type(const Point_3& query)
	{
		Vertex_handle v = nearest_boundary_vertex(query);
		if (v == Vertex_handle()) return false;

		std::cout << "vertex toggled " << v->point() << std::endl;
		toggle_vertex_type(v);

		return true;
	}

	Vertex_handle nearest_boundary_vertex(const Point_3& query)
	{
		double min_distance = 1e38;
        Point_2 query2d(query.x(), query.y());
		Vertex_handle nearest = Vertex_handle(); // aka NULL
		for (Vertex_iterator v = this->vertices_begin(); v != this->vertices_end(); v++)
		{
            Point_2 p = v->get_2d_point();
            const FT d2 = CGAL::squared_distance(p, query2d);
            if (d2 < min_distance)
            {
                min_distance = d2;
                nearest = v;
            }
        }
		return nearest;
	}

	FT sqd(const Point_3& a,
		const Point_3& b)
	{
		return CGAL::squared_distance(a,b);
	}

	FT distance(Vertex_handle v1,
		Vertex_handle v2)
	{
		const Point_3& a = v1->point();
		const Point_3& b = v2->point();
		return std::sqrt(sqd(a,b));
	}

	int index_boundary_vertices_per_boundary()
    {
        const int FREE = -1;
        const int DONE = +1;
        this->tag_halfedges(FREE);

		// init 
		this->set_all_vertex_boundary_indices(-1);
		
        int index = 0;
        for (Halfedge_iterator
             he  = this->halfedges_begin();
             he != this->halfedges_end();
             he++)
        {
            if (!he->is_border()) continue;
            if (he->tag() == DONE) continue;
            
            Halfedge_handle curr = he;
            do
            {
                curr->vertex()->boundary_index() = index;
                curr->tag() = DONE;
                curr = curr->next();
            }
            while(curr != he);
            index++;
        }
        return index;
    }

	// return index of longest boundary
	int index_longest_boundary()
	{
        const int FREE = -1;
        const int DONE = +1;
        this->tag_halfedges(FREE);

		FT max_len = 0.0;
        int index_longest = 0;
        int index = 0;
		double len = 0.0;
        for (Halfedge_iterator
             he  = this->halfedges_begin();
             he != this->halfedges_end();
             he++)
        {
            if (!he->is_border()) continue;
            if (he->tag() == DONE) continue;
            
            Halfedge_handle curr = he;
            do
            {
				len += len_3d(curr);
                curr->tag() = DONE;
                curr = curr->next();
            }
            while(curr != he);

			if(len > max_len)
			{
				index_longest = index;
				max_len = len;
			}

            index++;
			len = 0.0; // reset len
        }
        return index_longest;
	}


	void compute_generators(const int nb_boundaries)
	{
		std::cout << "compute generators" << std::endl;

		// we assume that all vertices have a proper boundary index

		// pick longest boundary as ref. one
		int ref_index = index_longest_boundary();

		m_generators.clear();
		for (int i=0; i<nb_boundaries; i++)
		{
			if (i != ref_index)
			{
				std::cout << "generators " << ref_index << " -> " << i << std::endl;
				Generator generator;
				bool ok = compute_generator(generator, ref_index, i, 10); // FERN: trial=10
				m_generators.push_back(generator);
                
                if (!ok)
                {
                    std::cout << "compute generator failed" << std::endl;
                }
			}
		}
	}

	// compute generator between boundary 1 and 2, and uses several random trials 
	// in order to retain the shortest one
	bool compute_generator(Generator& generator,
                           const int index1,
                           const int index2,
                           const int trial)
	{
		std::vector<Vertex_handle> vertices1;
		std::vector<Vertex_handle> vertices2;
		find_boundary_vertices_with_min_degree(vertices1, index1, 3);
		find_boundary_vertices_with_min_degree(vertices2, index2, 3);

		if (vertices1.empty() || vertices2.empty())
		{
			std::cerr << red << "unable to compute generator" << white << std::endl;
			return false;
		}

		srand(0);
		bool success = false;
		FT min_distance = 1e38;
		for (int i = 0; i < trial; i++)
		{
			Vertex_handle source = random_vertex(vertices1);
			Vertex_handle target = random_vertex(vertices2);
			FT distance = 0.0;
			//std::cout << "dijsktra " << source->index() << " -> " << target->index() << "...";
			bool ok = dijkstra_inner(source, target, distance);
			if (ok && distance < min_distance)
			{
				min_distance = distance;
				get_generator(target, generator); // this erases container generator
				success = true;
			}
			//std::cout << "done (distance " << distance << ")" << std::endl;
		}
		return success;
	}

	Vertex_handle random_vertex(std::vector<Vertex_handle>& vertices)
	{
		int index = random_int(0, vertices.size() - 1);
		return vertices[index];
	}

	// return arbitrary vertex with boundary index and of degree >= to min_degree
	void find_boundary_vertices_with_min_degree(std::vector<Vertex_handle>& vertices,
                                                const int boundary_index,
                                                const int min_degree)
	{
        for (Vertex_iterator v  = this->vertices_begin(); v != this->vertices_end(); v++)
        {
			if (is_border_vertex(v))
            {
				if (v->boundary_index() == boundary_index)
                {
					if ((int)v->vertex_degree() >= min_degree)
						vertices.push_back(v);
                }
            }
        }
	}

	// performs Dijstrka between two boundary vertices and
	// goes only through inner edges
	// returns true if succeeds and measure total distance
	// along path
	bool dijkstra_inner(Vertex_handle source,
                        Vertex_handle target,
                        FT& total_distance)
	{
		typedef std::priority_queue<Vertex_handle,
        std::vector<Vertex_handle>,
        Bigger_distance<Vertex_handle> > PQueue;

		if(!this->is_border_vertex(source) || !this->is_border_vertex(target))
		{
			std::cerr << red << "dijkstra_inner must run with two border vertices" 
				<< white << std::endl;
			return false;
		}

		if(source->vertex_degree() < 3 || target->vertex_degree() < 3)
		{
			std::cerr << red << "dijkstra_inner must run with two degree >= 3 vertices" 
				<< white << std::endl;
			return false;
		}

		PQueue queue;
		const int FREE = 0;
		const int DONE = 1;
		init_dijkstra(Vertex_handle(), FREE);

		// init priority queue
		source->distance() = 0.0f;
		queue.push(source);

		bool success = false;
		while(!queue.empty())
		{
			Vertex_handle v = queue.top();
			queue.pop();

			Halfedge_around_vertex_circulator he = v->vertex_begin();
			Halfedge_around_vertex_circulator end = he;
			CGAL_For_all(he,end)
			{
				// skip boundary edges
				if(is_border_edge(he)) continue;

				Vertex_handle nv = he->opposite()->vertex();
				double d = v->distance() + len_2d(he); // FERN: 2d or 3d?
				if (d < nv->distance())
				{
					nv->distance() = d;
					nv->from() = v;
				}
				if (v->tag() == FREE) queue.push(nv);
			}
            
			v->tag() = DONE;
			if (v == target)
			{
				success = true;
				break;
			}
		}

		// measure total distance along path
		if (success)
			total_distance = len_dijkstra(target);
		return success;
	}

	FT len_dijkstra(Vertex_handle target)
	{
		FT len = 0.0;
		Vertex_handle v = target;
		while(v->from() != Vertex_handle())
		{
			len += distance(v, v->from());
			v = v->from();
		}
		return len;
	}

    // from target to source
	void get_generator(Vertex_handle target, Generator& generator)
	{
		generator.clear();
		Vertex_handle v = target;
		while (v->from() != Vertex_handle())
		{
			Halfedge_handle he = get_halfedge_between(v, v->from());
			if (he != Halfedge_handle()) generator.push_back(he);
			v = v->from();
		}
	}

	void draw_generators()
	{
		typename std::vector<Generator>::iterator it;
		for(it = m_generators.begin(); it != m_generators.end(); it++)
			draw_generator(*it);
	}

	void draw_generator(Generator& gen)
	{
		::glBegin(GL_LINES);
		typename Generator::iterator it;
		for(it = gen.begin(); it != gen.end(); it++)
		{
			Halfedge_handle he = *it;
			const Point_3& p = get_source_vertex(he)->point();	
			const Point_3& q = get_target_vertex(he)->point();	
			::glVertex3d(p.x(), p.y(), 0.0); // p.z());
			::glVertex3d(q.x(), q.y(), 0.0); // q.z());
		}
		::glEnd();
	}

	// returns halfedge that goes from v1 to v2
	Halfedge_handle get_halfedge_between(Vertex_handle v1,  Vertex_handle v2)
	{
			Halfedge_around_vertex_circulator he = v1->vertex_begin();
			Halfedge_around_vertex_circulator end = he;
			CGAL_For_all(he, end)
			{
				Halfedge_handle ho = he->opposite();
				Vertex_handle nv = ho->vertex();
				if (nv == v2) return ho;
			}
			return Halfedge_handle();
	}

	void dijkstra(Vertex_handle source,
                  Vertex_handle target)
	{

		typedef typename std::priority_queue<Vertex_handle,
			std::vector<Vertex_handle>,
			Bigger_distance<Vertex_handle> > PQueue;
		PQueue queue;
		const int FREE = 0;
		const int DONE = 1;
		init_dijkstra(Vertex_handle(), FREE);

		// init priority queue
		source->distance() = 0.0f;
		queue.push(source);

		while(!queue.empty())
		{
			Vertex_handle v = queue.top();
			queue.pop();

			Halfedge_around_vertex_circulator he = v->vertex_begin();
			Halfedge_around_vertex_circulator end = he;
			CGAL_For_all(he,end)
			{
				Vertex_handle nv = he->opposite()->vertex();
				double d = v->distance() + len_3d(he);
				if(d < nv->distance())
				{
					nv->distance() = d;
					nv->from() = v;
				}
				if(v->tag() == FREE)
					queue.push(nv);
			}
			v->tag() = DONE;
			if(v == target)
				break; 
		}
	}

	void init_dijkstra(Vertex_handle root, const int init)
	{
		for(Vertex_iterator
            v = this->vertices_begin();
			v != this->vertices_end();
			v++)
		{
			v->distance() = 1e38f;
			v->from() = root;
			v->tag() = init;
		}
	}

	void gl_draw_dijkstra(Vertex_handle source,
		Vertex_handle target)
	{
		Vertex_handle v = target;
		::glColor3d(255,0,0);
		::glLineWidth(3.0f);
		::glBegin(GL_LINES);
		while (v->from() != Vertex_handle())
		{
			const Point_3& a = v->point();
			v = v->from();
			const Point_3& b = v->point();
			::glVertex3d(a.x(),a.y(),a.z());
			::glVertex3d(b.x(),b.y(),b.z());
		}
		::glEnd();
		::glFlush();
	}

	// normals (per	facet, then	per	vertex)
	void compute_normals()
	{
		compute_normals_per_facet();
		compute_normals_per_vertex();
	}
    
	void compute_normals_per_facet()
	{
		std::for_each(this->facets_begin(),this->facets_end(),Facet_normal());
	}
	
    void compute_normals_per_vertex()
	{
		std::for_each(this->vertices_begin(),this->vertices_end(),Vertex_normal());
	}
};

#endif
