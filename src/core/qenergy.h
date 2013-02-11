#ifndef _QENERGY_H_
#define _QENERGY_H_

#include "../third/energy.h"

template <class CENERGY>
class CQEnergy : public TemplateEnergy
{
public:
    typedef CENERGY Energy;
    typedef typename Energy::Mesh Mesh;
    
    typedef typename Mesh::FT         FT;
    typedef typename Mesh::Point_2    Point_2;
    typedef typename Mesh::Vector_2   Vector_2;
    typedef typename Mesh::Segment_2  Segment_2;
    typedef typename Mesh::Triangle_2 Triangle_2;
    typedef typename Mesh::Point_3    Point_3;
    typedef typename Mesh::Vector_3   Vector_3;
    typedef typename Mesh::Segment_3  Segment_3;
    typedef typename Mesh::Triangle_3 Triangle_3;

    typedef typename Mesh::Kernel Kernel;
    typedef typename Kernel::Tetrahedron_3 Tetrahedron_3;
    
    typedef typename Mesh::Edge_iterator Edge_iterator;
	typedef typename Mesh::Facet_handle   Facet_handle;
	typedef typename Mesh::Facet_iterator Facet_iterator;
	typedef typename Mesh::Vertex_handle   Vertex_handle;
	typedef typename Mesh::Vertex_iterator Vertex_iterator;
	typedef typename Mesh::Halfedge_handle   Halfedge_handle;
	typedef typename Mesh::Halfedge_iterator Halfedge_iterator;
	typedef typename Mesh::Halfedge_around_facet_circulator  Halfedge_around_facet_circulator;
	typedef typename Mesh::Halfedge_around_vertex_circulator Halfedge_around_vertex_circulator;
    
private:
    bool m_2d_mode;
    Energy* m_energy;
    Mesh* m_backup_mesh;
    
public:
    CQEnergy(Energy* energy)
    {
        m_2d_mode = true;
        m_energy = energy;
        m_backup_mesh = NULL;
    }
    
    ~CQEnergy()
    {
        if (m_backup_mesh) delete(m_backup_mesh);
    }
    
    bool get_2d_mode() const { return m_2d_mode; }
    void set_2d_mode(bool yes) { m_2d_mode = yes; }
    
    void set_backup()
    {
        std::cout << "set_backup...";
        if (m_backup_mesh) delete(m_backup_mesh);
        Mesh* mesh = m_energy->get_mesh();
        m_backup_mesh = m_energy->duplicate(mesh);
        std::cout << "done" << std::endl;
    }

    bool is_mesh_intrinsic_delaunay() const
    {
        Mesh* mesh = m_energy->get_mesh();
        for(Edge_iterator
			he  = mesh->edges_begin();
			he != mesh->edges_end();
            he++)
		{
            if (mesh->is_border_edge(he)) continue;
            if (mesh->get_surf_star1(he) < 0.0) return false;
        }
        return true;
    }
    
    double evaluate() const
    {
        Mesh* mesh = m_energy->get_mesh();
        FT sum = 0.0;
        for (Facet_iterator
             f  = mesh->facets_begin();
             f != mesh->facets_end();
             ++f)
        {
            sum += compute_cvt(f);
            //sum += compute_vol(f);
        }
        //sum += compute_inertia();
        return sum;
    }
    
    double compute_inertia() const
    {
        double sum = 0.0;
        Mesh* mesh = m_energy->get_mesh();
        for (Vertex_iterator
             v = mesh->vertices_begin();
             v != mesh->vertices_end();
             v++)
        {
            Vector_3 vec = v->get_3d_point() - v->backup();
            //sum += v->inertia()*vec.z()*vec.z();
            sum += v->inertia() * vec * vec;
        }
        return 0.5*m_energy->get_qinertia()*sum;
    }
    
    double compute_cvt(Facet_handle f) const
    {
        double sum = 0.0;
        Mesh* mesh = m_energy->get_mesh();
        Halfedge_around_facet_circulator he = f->facet_begin();
        Halfedge_around_facet_circulator end = he;
        CGAL_For_all(he, end)
        {
            FT len = mesh->len_3d(he);
            FT d = 0.5*len;
            FT h = len*0.5*mesh->cotangent_3d_from_h(he);
            // twice because dij == dji
            sum += std::pow(d,3)*h/2.0 + std::pow(h,3)*d/6.0;
        }
        return mesh->nb_vertices()*f->rho()*sum;
    }
    
    FT compute_vol(Facet_handle f) const
    {
        Mesh* mesh = m_energy->get_mesh();
        Triangle_3 t = mesh->triangle_3d(f);
        Tetrahedron_3 tet(t[0], t[2], t[1], CGAL::ORIGIN);
        return tet.volume();
    }
    
    Array gradient() const
    {
        Mesh* mesh = m_energy->get_mesh();
        Array G(3*mesh->nb_free_vertices());
        for (Vertex_iterator
             v = mesh->vertices_begin();
             v != mesh->vertices_end();
             v++)
        {
            if (m_energy->disregard_vertex(v)) continue;

            Vector_3 grad = compute_cvt_grad(v);
            //grad = grad + compute_vol_grad(v);
            //grad = grad + compute_inertia_grad(v);
            
            G[3*v->free_index()  ] = grad.x();
            G[3*v->free_index()+1] = grad.y();
            if (m_2d_mode) G[3*v->free_index()+2] = 0.0;
            else G[3*v->free_index()+2] = grad.z();
        }
        return G;
    }
    
    Vector_3 compute_inertia_grad(Vertex_handle v) const
    {
        Vector_3 vec = v->get_3d_point() - v->backup();
        //Vector_3 g(0.0, 0.0, vec.z());
        Vector_3 g = vec;
        return m_energy->get_qinertia() * v->inertia() * g;
    }
    
    Vector_3 compute_cvt_grad(Vertex_handle v) const
    {
        FT area = 0.0;
        Vector_3 bary(0.0,0.0,0.0);
        Vector_3 bdry(0.0,0.0,0.0);
        Mesh* mesh = m_energy->get_mesh();
        Halfedge_around_vertex_circulator hki = v->vertex_begin();
        Halfedge_around_vertex_circulator end = hki;
        CGAL_For_all(hki, end)
        {
            if (hki->is_border()) continue;
            
            Facet_handle f = hki->facet();
            Vector_3 n = mesh->get_3d_normal(f);
            
            const Point_3 pi = mesh->get_target_vertex(hki)->get_3d_point();
            const Point_3 pj = mesh->get_opposite_vertex(hki)->get_3d_point();
            const Point_3 pk = mesh->get_source_vertex(hki)->get_3d_point();
            
            Point_3 c = CGAL::circumcenter(pi, pj, pk);
            Point_3 pij = CGAL::midpoint(pi, pj);
            Point_3 pik = CGAL::midpoint(pi, pk);
            
            Triangle_3 t0(c, pik, pi);
            Triangle_3 t1(c, pi, pij);
            
            FT a0 = std::sqrt(t0.squared_area());
            FT a1 = std::sqrt(t1.squared_area());
            
            Vector_3 n0 = t0.supporting_plane().orthogonal_vector();
            Vector_3 n1 = t1.supporting_plane().orthogonal_vector();
            
            if (n*n0 < 0.0) a0 = -a0;
            if (n*n1 < 0.0) a1 = -a1;
                        
            FT lki = mesh->len_3d(hki);
            FT lij = mesh->len_3d(hki->next());
            
            Vector_3 nki = - mesh->get_3d_edge_normal(hki);
            Vector_3 nij = - mesh->get_3d_edge_normal(hki->next());

            Point_3 b0 = CGAL::centroid(t0);
            Point_3 b1 = CGAL::centroid(t1);
            
            area += f->rho()*( a0 + a1 );
            bdry = bdry + f->rho()*( std::pow(lij,3)*nij + std::pow(lki,3)*nki );
            bary = bary + f->rho()*( a0*(b0-CGAL::ORIGIN) + a1*(b1-CGAL::ORIGIN) );
        }
        Vector_3 p = v->get_3d_point() - CGAL::ORIGIN;
        return mesh->nb_vertices()*( 2.0*(area*p - bary) + (1.0/24.0)*bdry );
    }

    Vector_3 compute_vol_grad(Vertex_handle v) const
    {
        Vector_3 sum_vec(0.0,0.0,0.0);
        Mesh* mesh = m_energy->get_mesh();
        Halfedge_around_vertex_circulator hki = v->vertex_begin();
        Halfedge_around_vertex_circulator end = hki;
        CGAL_For_all(hki, end)
        {
            if (hki->is_border()) continue;
            Vector_3 pj = mesh->get_opposite_vertex(hki)->get_3d_point() - CGAL::ORIGIN;
            Vector_3 pk = mesh->get_source_vertex(hki)->get_3d_point() - CGAL::ORIGIN;
            sum_vec = sum_vec + CGAL::cross_product(pj, pk)/6.0;
        }
        return sum_vec;
    }

    Array get_variables() const
    {
        Mesh* mesh = m_energy->get_mesh();
        Array X(3*mesh->nb_free_vertices());
        for (Vertex_iterator
             v = mesh->vertices_begin();
             v != mesh->vertices_end();
             v++)
        {
            if (m_energy->disregard_vertex(v)) continue;
            Point_3 p = v->get_3d_point();
            X[3*v->free_index()  ] = p.x();
            X[3*v->free_index()+1] = p.y();
            X[3*v->free_index()+2] = p.z();
        }
        return X;
    }
    
    bool set_variables(const Array& xyz)
    {
        std::vector<Point_3> p = convert_array_to_vector(xyz);
        if (m_energy->get_fixed_connectivity()) return move_mesh(p);
        return move_and_adapt_mesh(p);
    }

    bool move_mesh(const std::vector<Point_3>& points)
    {
        Mesh* mesh = m_energy->get_mesh();
        for (Vertex_iterator
             v = mesh->vertices_begin();
             v != mesh->vertices_end();
             v++)
        {
            if (m_energy->disregard_vertex(v)) continue;
            const Point_3 p = points[v->free_index()];
            const Point_2 q(p.x(), p.y());
            if (mesh->valid_embedding(v, q))
                v->set_3d_point(p);
        }
        m_energy->pre_compute_data();
        if (!m_energy->get_validity()) return true;
        return is_mesh_intrinsic_delaunay();
    }

    bool move_and_adapt_mesh(const std::vector<Point_3>& points)
    {
        m_energy->duplicate_and_set(m_backup_mesh);
        Mesh* mesh = m_energy->get_mesh();
        Meshing meshing(mesh);
        for (Vertex_iterator
             v = mesh->vertices_begin();
             v != mesh->vertices_end();
             v++)
        {
            if (m_energy->disregard_vertex(v)) continue;
            const Point_3 p = points[v->free_index()];
            const Point_2 q(p.x(), p.y());
            if (mesh->valid_embedding(v, q))
                v->set_3d_point(p);
            meshing.stack_intrinsic_delaunay_flip(v);
        }
        meshing.intrinsic_delaunay_flip_edges();
        m_energy->pre_compute_data();
        if (!m_energy->get_validity()) return true;
        return is_mesh_intrinsic_delaunay();
    }
    
    Array convert_vector_to_array(const std::vector<Point_3>& q) const
    {
        Array xyz(3*q.size());
        for (unsigned i = 0; i < q.size(); ++i)
        {
            const Point_3& p = q[i];
            xyz[3*i  ] = p.x();
            xyz[3*i+1] = p.y();
            xyz[3*i+2] = p.z();
        }
        return xyz;
    }
    
    std::vector<Point_3> convert_array_to_vector(const Array& xyz) const
    {
        std::vector<Point_3> q(xyz.size()/3);
        for (unsigned i = 0; i < q.size(); i++)
            q[i] = Point_3(xyz[3*i], xyz[3*i+1], xyz[3*i+2]);
        return q;
    }
    
    // debug //
    
    void debug_gradient()
    {
        Mesh* mesh = m_energy->get_mesh();
        FT epsilon = 1.0e-6;
        for (Vertex_iterator
             v  = mesh->vertices_begin();
             v != mesh->vertices_end();
             ++v)
        {
            if (mesh->is_vertex_grounded(v)) continue;
            Vector_3 grad = compute_cvt_grad(v) + compute_vol_grad(v);
            Vector_3 num_grad = compute_numerical_gradient(v, epsilon);
            Vector_3 res_vec = num_grad - grad;
            FT res = sqrt(res_vec*res_vec);
            if (res > epsilon)
            {
                std::cout << "Vertex" << v->index() << ": "
                << " res = " << res
                << " (g = " << grad << " ; num = " << num_grad << ")" << std::endl;
            }
        }
    }
    
    Vector_3 compute_numerical_gradient(Vertex_handle vi, FT epsilon)
    {
        FT gx = 0.0;
        FT gy = 0.0;
        FT gz = 0.0;
        Vector_3 dx(epsilon, 0.0, 0.0);
        Vector_3 dy(0.0, epsilon, 0.0);
        Vector_3 dz(0.0, 0.0, epsilon);
        
        Mesh* mesh = m_energy->get_mesh();
        Point_3 pi = vi->get_3d_point();
        for (unsigned i = 0; i < 3; ++i)
        {
            if (i == 0) vi->set_3d_point(pi + dx);
            if (i == 1) vi->set_3d_point(pi + dy);
            if (i == 2) vi->set_3d_point(pi + dz);
            FT pos = evaluate();
            
            if (i == 0) vi->set_3d_point(pi - dx);
            if (i == 1) vi->set_3d_point(pi - dy);
            if (i == 2) vi->set_3d_point(pi - dz);
            FT neg = evaluate();
            
            vi->set_3d_point(pi);
            FT grad = (pos - neg)/(2*epsilon);
            
            if (i == 0) gx = grad;
            if (i == 1) gy = grad;
            if (i == 2) gz = grad;
        }
        return Vector_3(gx, gy, gz);
    }
};

#endif
