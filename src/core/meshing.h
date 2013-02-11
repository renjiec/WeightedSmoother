#ifndef _MESHING_H_
#define _MESHING_H_

#include <set>
#include "phalfedge.h"

template <class CMesh>
class CMeshing
{
public:
    typedef CMesh Mesh;
    
	typedef typename Mesh::Kernel Kernel;
	typedef typename Kernel::FT   FT;
	typedef typename Kernel::Point_2 Point_2;
	typedef typename Kernel::Point_3 Point_3;
    
	typedef typename Mesh::Facet_handle    Facet_handle;
	typedef typename Mesh::Vertex_handle   Vertex_handle;
	typedef typename Mesh::Vertex_iterator Vertex_iterator;
	typedef typename Mesh::Edge_iterator   Edge_iterator;
	typedef typename Mesh::Halfedge_handle Halfedge_handle;
	typedef typename Mesh::Halfedge_around_vertex_circulator HV_circulator;
    
	typedef CPHalfedge<FT, Halfedge_handle> PHalfedge;
	typedef CDPQueue_short<PHalfedge> DPQueue_short;
	typedef CDPQueue_long<PHalfedge> DPQueue_long;
    
private:
	Mesh* m_mesh;
    
public:
	CMeshing(Mesh* pRemesh)
	{
		m_mesh = pRemesh;
	}
    
	~CMeshing() { }

	

    
    // split //
    
	template <class Feature_tree, class Oracle>
	unsigned split_long_edges(const FT& max_len,
                              Feature_tree& feature_tree,
							  Oracle* pOracle,
                              bool lift_on = true)
	{
		DPQueue_long long_edges;
		fill_queue_with_long_edges(long_edges, max_len);
        
		unsigned nb_split = 0;
		while (!long_edges.empty())
		{
			PHalfedge edge = long_edges.top();
			long_edges.pop();
            Halfedge_handle he = edge.halfedge();
            
            long_edges.remove(PHalfedge(he));
            long_edges.remove(PHalfedge(he->opposite()));
            
            Halfedge_handle hnew = split_edge(he, feature_tree, pOracle, lift_on);
            
            add_circular_long_edges(hnew->vertex(), long_edges, max_len);
			nb_split++;
        }
        return nb_split;
    }
    
	template <class Feature_tree, class Oracle>
    Halfedge_handle split_edge(Halfedge_handle he,
                               Feature_tree& feature_tree,
                               Oracle* pOracle,
                               bool lift_on)
    {
        // memorize source and target point
        const Point_3& source = m_mesh->get_source_vertex(he)->point();
        const Point_3& target = m_mesh->get_target_vertex(he)->point();
        const FT init_len = distance_3d(source, target);
        
        // set split point to midpoint by default
        const Point_3 midpoint_3d = m_mesh->midpoint_3d(he);
        Point_3 new_point = midpoint_3d;
        
        if (lift_on)
        {
            const FT z = pOracle->evaluate_height(to_2d(midpoint_3d));
            new_point = Point_3(midpoint_3d.x(), midpoint_3d.y(), z);
        }
        
        bool boundary = he->is_border() || he->opposite()->is_border();
        Halfedge_handle hnew = m_mesh->split_edge(he);
        hnew->vertex()->point() = new_point;
        
        // check that edge has been shortened to avoid infinite loop
        // and check valid embedding
        if (boundary) new_point = feature_tree.closest_point(new_point);
        if (new_point != source && new_point != target)
        {
            const FT len1 = distance_3d(new_point, source);
            const FT len2 = distance_3d(new_point, target);
            if (len1 < init_len && len2 < init_len)
                if (m_mesh->valid_embedding(hnew->vertex(), to_2d(new_point)))
                    hnew->vertex()->point() = new_point;
        }
        
        // split incident facets
        if (!hnew->is_border())
            m_mesh->split_facet(hnew, hnew->next()->next());
        
        if (!hnew->opposite()->is_border())
            m_mesh->split_facet(hnew->opposite()->next(),
                                hnew->opposite()->next()->next()->next());
        
        return hnew;
	}
    
	void fill_queue_with_long_edges(DPQueue_long& queue, const FT max_len)
	{
		for(Edge_iterator
			he = m_mesh->edges_begin();
			he != m_mesh->edges_end();
            he++)
		{
			const double len = m_mesh->len_3d(he);
			if (len > max_len) queue.push(PHalfedge(he, len));
		}
	}
    
	void add_circular_long_edges(Vertex_handle v, DPQueue_long& queue, const FT max_len)
	{
		HV_circulator he = v->vertex_begin();
		HV_circulator end = he;
		CGAL_For_all(he, end)
		{
			const double len = m_mesh->len_3d(he);
			if (len > max_len) queue.push(PHalfedge(he, len));
		}
	}
    
    // refine //
    
	template <class Facet_tree, class Feature_tree, class Oracle>
	unsigned refine_for_height(Facet_tree& facet_tree,
                               Feature_tree& feature_tree,
                               Oracle *pOracle,
                               const unsigned nb_probes,
                               const FT max_height)
	{
		DPQueue_long long_edges;
		fill_queue_with_long_height_edges(long_edges, pOracle, nb_probes, max_height);
        
		unsigned nb_split = 0;
		while (!long_edges.empty())
		{
			PHalfedge edge = long_edges.top();
			long_edges.pop();
            
			Halfedge_handle he = edge.halfedge();
			bool boundary = he->is_border() || he->opposite()->is_border();
            
			long_edges.remove(PHalfedge(he));
			long_edges.remove(PHalfedge(he->opposite()));
            
			// compute split point to midpoint by default
			const Point_2 midpoint = m_mesh->midpoint_2d(he);
			const FT z = pOracle->evaluate_height(midpoint);
			Point_3 lifted_midpoint(midpoint.x(), midpoint.y(), z);
            
			// special case: relocate to boundary
			bool do_split = !boundary;
            
			if(boundary)
			{
				// nearest boundary edges
				Point_3 new_pos = feature_tree.closest_point(lifted_midpoint);
				Point_2 new_pos2d = to_2d(new_pos);
                
				const Point_3& p = he->vertex()->point();
				const Point_3& q = he->opposite()->vertex()->point();
                
				// check that new position is different from edge end points
				if(new_pos != p && new_pos != q)
                    // check valid embedding and lift
                    if (m_mesh->valid_embedding_upon_boundary_split(he, new_pos2d))
                    {
                        const FT z = pOracle->evaluate_height(new_pos2d);
                        lifted_midpoint = Point_3(new_pos.x(), new_pos.y(), z);
                        do_split = true;
                    }
			}
            
			if(!do_split)
				continue;
            
			// perform split
			Halfedge_handle hnew = m_mesh->split_edge(he);
			hnew->vertex()->point() = lifted_midpoint;
			nb_split++;
            
			// split incident facets
			if (!hnew->is_border())
				m_mesh->split_facet(hnew, hnew->next()->next());
            
			if (!hnew->opposite()->is_border())
				m_mesh->split_facet(hnew->opposite()->next(),
                                    hnew->opposite()->next()->next()->next());
            
			// update queue
			add_circular_long_height_edges(hnew->vertex(),
                                           long_edges, pOracle, nb_probes, max_height);
            
            // debug
			if(nb_split > 1e4)
			{
				std::cerr << red << " premature loop ending in refine_for_height" << white << std::endl;
				break;
			}
            //
		}
		return nb_split;
	}
    
	void fill_queue_with_long_height_edges(DPQueue_long& queue,
                                           Oracle *pOracle,
                                           const unsigned nb_probes,
                                           const FT max_height)
	{
		for (Edge_iterator he = m_mesh->edges_begin(); he != m_mesh->edges_end(); he++)
        {
            FT error = m_mesh->evaluate_height_error(he, pOracle, nb_probes);
			if (error > max_height)
				queue.push(PHalfedge(he, error));
			else
				if(m_mesh->boundary_edge_to_conform(he))
					queue.push(PHalfedge(he, 1e6));
        }
	}
    
	void add_circular_long_height_edges(Vertex_handle v, DPQueue_long& queue,
                                        Oracle *pOracle,
                                        const unsigned nb_probes,
                                        const FT max_height)
	{
		HV_circulator he = v->vertex_begin();
		HV_circulator end = he;
		CGAL_For_all(he, end)
        {
            FT error = m_mesh->evaluate_height_error(he, pOracle, nb_probes);
			if (error > max_height)
				queue.push(PHalfedge(he, error));
			else
				if(m_mesh->boundary_edge_to_conform(he))
					queue.push(PHalfedge(he, 1e6));
        }
	}
    
    // collapse //
    
	template <class Feature_tree, class Oracle>
	unsigned collapse_short_edges(const FT min_len,
                                  Feature_tree& feature_tree,
								  Oracle* pOracle)
	{
		DPQueue_short short_edges;
		fill_queue_with_collapsible_short_edges(short_edges, min_len);
        
		unsigned nb_collapse = 0;
		while (!short_edges.empty())
		{
			PHalfedge edge = short_edges.top();
			short_edges.pop();
            
			Halfedge_handle he = edge.halfedge();
			if (!m_mesh->is_collapsible(he))
			{
				if(he->opposite()->is_border())
					std::cerr << "non-collapsible boundary edge popped: " << &*he << std::endl;
				else
					std::cerr << "non-collapsible inner edge popped: " << &*he << std::endl;
				continue;
			}

			// get vertices
			Vertex_handle vq = m_mesh->get_target_vertex(he);
			Vertex_handle vp = m_mesh->get_source_vertex(he);

			// compute point location
			bool do_collapse = false;
			bool do_tag = false;
			Point_3 final_location;
			if(m_mesh->is_vertex_grounded(vp) || m_mesh->is_vertex_grounded(vq))
			{
				// p is not grounded and q is grounded -> q
				if(!m_mesh->is_vertex_grounded(vp) && m_mesh->is_vertex_grounded(vq))
				{
					const Point_3& q = vq->point();
					if (m_mesh->valid_embedding_upon_collapse(he, to_2d(q)))
					{
						final_location = q;
						do_collapse = true;
						do_tag = true;
					}
				}
				// p is grounded and q is not grounded -> p
				if(m_mesh->is_vertex_grounded(vp) && !m_mesh->is_vertex_grounded(vq))
				{
					const Point_3& p = vp->point();
					if (m_mesh->valid_embedding_upon_collapse(he, to_2d(p)))
					{
						final_location = p;
						do_collapse = true;
						do_tag = true;
					}
				}
			}
			else
			{
				// compute new point location
				Point_3 midpoint_3d = m_mesh->midpoint_3d(he);
				const Point_2 midpoint_2d = to_2d(midpoint_3d);

				const FT z = pOracle->evaluate_height(midpoint_2d);
				Point_3 lifted_midpoint(midpoint_2d.x(), midpoint_2d.y(), z);

				final_location = lifted_midpoint;
				if (m_mesh->touch_border(he))
					final_location = feature_tree.closest_point(midpoint_3d);
				//else
				//	final_location = facet_tree.closest_point(midpoint3d);
				Point_2 projected_midpoint2d = to_2d(final_location);

				// check valid embedding with projected_midpoint2d
				if (m_mesh->valid_embedding_upon_collapse(he, projected_midpoint2d))
					do_collapse = true;
			}
            
			if(!do_collapse)
				continue;
            
			// remove edges from p. queue before collapse
			remove_link_short_edges(vp, short_edges);
			remove_link_short_edges(vq, short_edges);
            
			// join two pairs of facets before collapse
 			if (!m_mesh->join_facets_before_collapse(he))
			{
				std::cerr << red << "unable to join faces (premature ending)" << white << std::endl;
				return nb_collapse;
			}
            
			// collapse
			Halfedge_handle he_joined = m_mesh->join_vertex(he);
			Vertex_handle v_joined = he_joined->vertex();
			assert(v_joined == vq);
			v_joined->point() = final_location;

			// in case of collapse with fixed tag
			if(do_tag)
				v_joined->type() = Mesh::FIXED;

			nb_collapse++;
            
			// add edges to queue
			add_link_collapsible_short_edges(v_joined, short_edges, min_len);
		}
        
		return nb_collapse;
	}
    
	void remove_link_short_edges(Vertex_handle v, DPQueue_short& queue)
	{
		HV_circulator he = v->vertex_begin();
		HV_circulator end = he;
		CGAL_For_all(he, end)
		{
			Vertex_handle v = m_mesh->get_source_vertex(he);
			remove_incident_short_edges(v, queue);
		}
	}
    
	void remove_incident_short_edges(Vertex_handle v, DPQueue_short& queue)
	{
		HV_circulator he = v->vertex_begin();
		HV_circulator end = he;
		CGAL_For_all(he, end)
		{
			queue.remove(PHalfedge(he));
			queue.remove(PHalfedge(he->opposite()));
		}
	}
    
	void fill_queue_with_collapsible_short_edges(DPQueue_short& queue, const FT min_len)
	{
		for(Edge_iterator
			he = m_mesh->edges_begin();
			he != m_mesh->edges_end();
            he++)
		{
			add_to_collapsible_short_edges(he, queue, min_len);
		}
	}
    
	void collect_link_short_inner_halfedges(Vertex_handle v, std::set<Halfedge_handle>& halfedges, const FT min_len)
	{
		HV_circulator he = v->vertex_begin();
		HV_circulator end = he;
		CGAL_For_all(he, end)
		{
			Vertex_handle v = m_mesh->get_source_vertex(he);
			collect_incident_short_inner_halfedges(v, halfedges, min_len);
		}
	}
	void collect_incident_short_inner_halfedges(Vertex_handle v,
                                                std::set<Halfedge_handle>& halfedges,
                                                const FT min_len)
	{
		HV_circulator he = v->vertex_begin();
		HV_circulator end = he;
		CGAL_For_all(he, end)
		{
			const double len = m_mesh->len_3d(he);
			if (len > min_len) continue;
            
			Halfedge_handle ho = he->opposite();
			Halfedge_handle to_add = he->is_border() ? ho : he;
            
			// add if not already added (nor its opposite)
			if(halfedges.find(ho) == halfedges.end())
				if(halfedges.find(he) == halfedges.end())
					halfedges.insert(to_add);
		}
	}
    
	void add_link_collapsible_short_edges(Vertex_handle v, DPQueue_short& queue, const FT min_len)
	{
		std::set<Halfedge_handle> halfedges;
		collect_link_short_inner_halfedges(v, halfedges, min_len);
        
		typename std::set<Halfedge_handle>::iterator it;
		for(it = halfedges.begin(); it != halfedges.end(); it++)
		{
			Halfedge_handle he = *it;
			const double len = m_mesh->len_3d(he);
			if (m_mesh->is_collapsible(he))
				queue.push(PHalfedge(he, len));
		}
	}
    
	// if shorter than min_len, add he (or its opposite if on boundary)
	// to given priority queue
	void add_to_collapsible_short_edges(Halfedge_handle he, DPQueue_short& queue, const FT min_len)
	{
		// measure edge len
		const double len = m_mesh->len_3d(he);
		if (len > min_len) return;
        
		// do not add border half edges
		if (!he->is_border())
		{
			if (m_mesh->is_collapsible(he))
				queue.push(PHalfedge(he, len));
		}
		else
		{
			Halfedge_handle ho = he->opposite();
			if (m_mesh->is_collapsible(ho))
				queue.push(PHalfedge(ho, len));
		}
		//if (m_mesh->is_vertex_grounded(h->opposite()->vertex())) continue;
	}
    
    // remove 3-valence vertices //
    
	unsigned remove_degree3_vertices()
	{
		unsigned removed = 0;
		Vertex_handle v = Vertex_handle();
		while(find_inner_degree3_vertex(v))
		{
			m_mesh->erase_center_vertex(v->halfedge());
			removed++;
		}
		return removed;
	}
    
	bool find_inner_degree3_vertex(Vertex_handle& vh)
	{
		for (Vertex_iterator
             v = m_mesh->vertices_begin();
             v != m_mesh->vertices_end();
             v++)
		{
			if (m_mesh->is_border_vertex(v)) continue;
			if (v->degree() != 3) continue;
            
            vh = v;
            return true;
		}
		return false;
	}
    
    // regular_flip //
    
    // flip edges while seeding only around v
	unsigned int stack_regular_flip(Vertex_handle v)
	{
		DPQueue_long edges;
		fill_queue_with_incident_regular_flippable_edges(v, edges);
        
		unsigned nb_flip = 0;
		while (!edges.empty())
		{
			PHalfedge edge = edges.top();
			edges.pop();
			Halfedge_handle he = edge.halfedge();
            
			remove_edges_before_flip(edges, he);
			Halfedge_handle hnew = m_mesh->flip_edge(he);
			nb_flip++;
			add_regular_flippable_edges_after_flip(edges, hnew);
		}
		return nb_flip;
	}
    
	unsigned regular_flip_edges()
	{
		DPQueue_long edges;
		fill_queue_with_regular_flippable_edges(edges);
        
		unsigned nb_flip = 0;
		while (!edges.empty())
		{
			PHalfedge edge = edges.top();
			edges.pop();
			Halfedge_handle he = edge.halfedge();
            
			remove_edges_before_flip(edges, he);
			Halfedge_handle hnew = m_mesh->flip_edge(he);
			nb_flip++;
			add_regular_flippable_edges_after_flip(edges, hnew);
            
            // HACK
            //break;
            //
		}
		return nb_flip;
	}
    
	void fill_queue_with_regular_flippable_edges(DPQueue_long& queue)
	{
		for(Edge_iterator
			he = m_mesh->edges_begin();
			he != m_mesh->edges_end();
            he++)
		{
			if (!m_mesh->is_regular_flippable(he)) continue;
			const double key = m_mesh->get_star1(he);
            //const double key = m_mesh->len_2d(he);
			queue.push(PHalfedge(he, key));
		}
	}
    
    void fill_queue_with_incident_regular_flippable_edges(Vertex_handle v,
                                                          DPQueue_long& queue)
	{
		HV_circulator he = v->vertex_begin();
		HV_circulator end = he;
		CGAL_For_all(he, end)
		{
			if (!m_mesh->is_regular_flippable(he)) continue;
			const double key = m_mesh->get_star1(he);
            //const double key = m_mesh->len_2d(he);
			queue.push(PHalfedge(he, key));
		}
	}
    
    void add_regular_flippable_edges_after_flip(DPQueue_long& queue, Halfedge_handle he)
	{
		add_link_regular_flippable_edges(he->vertex(), queue);
		add_link_regular_flippable_edges(he->next()->vertex(), queue);
		add_link_regular_flippable_edges(he->opposite()->vertex(), queue);
		add_link_regular_flippable_edges(he->opposite()->next()->vertex(), queue);
	}
    
	void add_link_regular_flippable_edges(Vertex_handle v, DPQueue_long& queue)
	{
		HV_circulator he = v->vertex_begin();
		HV_circulator end = he;
		CGAL_For_all(he, end)
		{
			add_regular_flippable_edge_after_flip(queue, he);         // incident edge
			add_regular_flippable_edge_after_flip(queue, he->prev()); // link boundary
		}
	}
    
	bool add_regular_flippable_edge_after_flip(DPQueue_long& queue, Halfedge_handle he)
	{
		// exit if queue contains he or opposite
		if (queue.contains(PHalfedge(he))) return false;
		if (queue.contains(PHalfedge(he->opposite()))) return false;
        
		if (!m_mesh->is_regular_flippable(he)) return false;
        const double key = m_mesh->get_star1(he);
        //const double key = m_mesh->len_2d(he);
        queue.push(PHalfedge(he, key));
		return true;
	}
    
    // delaunay_flip //
    
	// flip edges while seeding only around v
	unsigned int stack_delaunay_flip(Vertex_handle v)
	{
		DPQueue_long edges;
		fill_queue_with_incident_delaunay_flippable_edges(v, edges);
        
		unsigned nb_flip = 0;
		while (!edges.empty())
		{
			PHalfedge edge = edges.top();
			edges.pop();
			Halfedge_handle he = edge.halfedge();
            
			remove_edges_before_flip(edges, he);
			Halfedge_handle hnew = m_mesh->flip_edge(he);
			nb_flip++;
			add_delaunay_flippable_edges_after_flip(edges, hnew);
		}
		return nb_flip;
	}
    
	// intrinsic_delaunay flip edges while seeding only around v
	unsigned int stack_intrinsic_delaunay_flip(Vertex_handle v)
	{
		DPQueue_long edges;
		fill_queue_with_incident_intrinsic_delaunay_flippable_edges(v, edges);
        
		unsigned nb_flip = 0;
		while (!edges.empty())
		{
			PHalfedge edge = edges.top();
			edges.pop();
			Halfedge_handle he = edge.halfedge();
            
			// FIXME
			if(nb_flip > 1e3)
			{
				std::cout << "premature ending" << std::endl;
				break;
			}
            
			remove_edges_before_flip(edges, he);
			Halfedge_handle hnew = m_mesh->flip_edge(he);
			nb_flip++;
			add_intrinsic_delaunay_flippable_edges_after_flip(edges, hnew);
		}
		return nb_flip;
	}
    
	unsigned intrinsic_delaunay_flip_edges()
	{
		DPQueue_long edges;
		fill_queue_with_intrinsic_delaunay_flippable_edges(edges);
        
		unsigned nb_flip = 0;
		while (!edges.empty())
		{
			PHalfedge edge = edges.top();
			edges.pop();
			Halfedge_handle he = edge.halfedge();
            
			// FIXME
			if(nb_flip > 1e3)
			{
				std::cout << "premature ending" << std::endl;
				break;
			}
            
			remove_edges_before_flip(edges, he);
			Halfedge_handle hnew = m_mesh->flip_edge(he);
			nb_flip++;
			add_intrinsic_delaunay_flippable_edges_after_flip(edges, hnew);
		}
		return nb_flip;
	}
    
    
	unsigned delaunay_flip_edges()
	{
		DPQueue_long edges;
		fill_queue_with_delaunay_flippable_edges(edges);
        
		unsigned nb_flip = 0;
		while (!edges.empty())
		{
			PHalfedge edge = edges.top();
			edges.pop();
			Halfedge_handle he = edge.halfedge();
            
			remove_edges_before_flip(edges, he);
			Halfedge_handle hnew = m_mesh->flip_edge(he);
			nb_flip++;
			add_delaunay_flippable_edges_after_flip(edges, hnew);
		}
		return nb_flip;
	}
    
	void fill_queue_with_intrinsic_delaunay_flippable_edges(DPQueue_long& queue)
	{
		for(Edge_iterator
			he = m_mesh->edges_begin();
			he != m_mesh->edges_end();
            he++)
		{
			if (!m_mesh->is_intrinsic_delaunay_flippable(he)) continue;
			const double key = m_mesh->angle_increase_upon_flip(he);
			queue.push(PHalfedge(he, key));
		}
	}
    
    
	void fill_queue_with_delaunay_flippable_edges(DPQueue_long& queue)
	{
		for(Edge_iterator
			he = m_mesh->edges_begin();
			he != m_mesh->edges_end();
            he++)
		{
			if (!m_mesh->is_delaunay_flippable(he)) continue;
			const double key = m_mesh->get_unweighted_star1(he);
			queue.push(PHalfedge(he, key));
		}
	}
    
	void fill_queue_with_incident_delaunay_flippable_edges(Vertex_handle v,
                                                           DPQueue_long& queue)
	{
		HV_circulator he = v->vertex_begin();
		HV_circulator end = he;
		CGAL_For_all(he, end)
		{
			if (!m_mesh->is_delaunay_flippable(he)) continue;
			const double key = m_mesh->get_unweighted_star1(he);
			queue.push(PHalfedge(he, key));
		}
	}
    
	void fill_queue_with_incident_intrinsic_delaunay_flippable_edges(Vertex_handle v,
                                                                     DPQueue_long& queue)
	{
		HV_circulator he = v->vertex_begin();
		HV_circulator end = he;
		CGAL_For_all(he, end)
		{
			if (!m_mesh->is_delaunay_flippable(he)) continue;
			const double key = m_mesh->angle_increase_upon_flip(he);
			queue.push(PHalfedge(he, key));
		}
	}
    
	void add_delaunay_flippable_edges_after_flip(DPQueue_long& queue, Halfedge_handle he)
	{
		add_link_delaunay_flippable_edges(he->vertex(), queue);
		add_link_delaunay_flippable_edges(he->next()->vertex(), queue);
		add_link_delaunay_flippable_edges(he->opposite()->vertex(), queue);
		add_link_delaunay_flippable_edges(he->opposite()->next()->vertex(), queue);
	}
    
	void add_link_delaunay_flippable_edges(Vertex_handle v, DPQueue_long& queue)
	{
		HV_circulator he = v->vertex_begin();
		HV_circulator end = he;
		CGAL_For_all(he, end)
		{
			add_delaunay_flippable_edge_after_flip(queue, he);         // incident edge
			add_delaunay_flippable_edge_after_flip(queue, he->prev()); // link boundary
		}
	}
    
	bool add_delaunay_flippable_edge_after_flip(DPQueue_long& queue, Halfedge_handle he)
	{
		// exit if queue contains he or opposite
		if (queue.contains(PHalfedge(he))) return false;
		if (queue.contains(PHalfedge(he->opposite()))) return false;
        
		if (!m_mesh->is_delaunay_flippable(he)) return false;
        const double key = m_mesh->get_unweighted_star1(he);
        queue.push(PHalfedge(he, key));
		return true;
	}
    
	void add_intrinsic_delaunay_flippable_edges_after_flip(DPQueue_long& queue, Halfedge_handle he)
	{
		add_link_intrinsic_delaunay_flippable_edges(he->vertex(), queue);
		add_link_intrinsic_delaunay_flippable_edges(he->next()->vertex(), queue);
		add_link_intrinsic_delaunay_flippable_edges(he->opposite()->vertex(), queue);
		add_link_intrinsic_delaunay_flippable_edges(he->opposite()->next()->vertex(), queue);
	}
    
	void add_link_intrinsic_delaunay_flippable_edges(Vertex_handle v, DPQueue_long& queue)
	{
		HV_circulator he = v->vertex_begin();
		HV_circulator end = he;
		CGAL_For_all(he, end)
		{
			add_intrinsic_delaunay_flippable_edge_after_flip(queue, he);         // incident edge
			add_intrinsic_delaunay_flippable_edge_after_flip(queue, he->prev()); // link boundary
		}
	}
    
	bool add_intrinsic_delaunay_flippable_edge_after_flip(DPQueue_long& queue, Halfedge_handle he)
	{
		// exit if queue contains he or opposite
		if (queue.contains(PHalfedge(he))) return false;
		if (queue.contains(PHalfedge(he->opposite()))) return false;
        
		if (!m_mesh->is_intrinsic_delaunay_flippable(he)) return false;
        const double key = m_mesh->angle_increase_upon_flip(he);
        queue.push(PHalfedge(he, key));
		return true;
	}
    
	void remove_edges_before_flip(DPQueue_long& queue, Halfedge_handle he)
	{
		remove_link_edges(he->vertex(), queue);
		remove_link_edges(he->next()->vertex(), queue);
		remove_link_edges(he->opposite()->vertex(), queue);
		remove_link_edges(he->opposite()->next()->vertex(), queue);
	}
    
	void remove_link_edges(Vertex_handle v, DPQueue_long& queue)
	{
		HV_circulator he = v->vertex_begin();
		HV_circulator end = he;
		CGAL_For_all(he, end)
		{
			queue.remove(PHalfedge(he));
			queue.remove(PHalfedge(he->opposite()));
            
			queue.remove(PHalfedge(he->prev()));
			queue.remove(PHalfedge(he->prev()->opposite()));
		}
	}
    
    // smooth //
    
	template <class Surface_tree, class Feature_tree, class Oracle>
	void smooth(Surface_tree& height_tree,
		        Feature_tree& feature_tree,
		        Oracle *pOracle,
                const int iter)
	{
		for(int i = 0; i < iter; i++)
			laplacian_smoothing(height_tree, feature_tree, pOracle);
	}
    
	template <class Surface_tree, class Feature_tree, class Oracle>
	void laplacian_smoothing(Surface_tree& height_tree,
                             Feature_tree& feature_tree, Oracle *pOracle)
	{
		for (Vertex_iterator
             v = m_mesh->vertices_begin();
             v != m_mesh->vertices_end();
             v++)
		{
			// skip fixed vertices
			if(m_mesh->is_vertex_grounded(v))
				continue;

			bool boundary = m_mesh->is_border_vertex(v);
			Point_3 centroid3d = m_mesh->barycenter_3d(v, boundary);
			Point_2 centroid2d = to_2d(centroid3d);
            
			Point_3 lifted_centroid = height_tree.closest_point(centroid3d);
            
			//const FT z = pOracle->evaluate_height(centroid2d);
			//Point_3 lifted_centroid(centroid2d.x(), centroid2d.y(), z);
            
			Point_3 new_point3d = lifted_centroid;
			if (boundary)
				new_point3d = feature_tree.closest_point(centroid3d);
			Point_2 new_point2d = to_2d(new_point3d);
            
			if (m_mesh->valid_embedding(v, new_point2d))
				v->point() = new_point3d;
		}
	}


	// ODT smoothing

	template <class Oracle>
	void odt_2d(Oracle* pOracle, const int iter)
	{
		for(int i = 0; i < iter; i++)
			one_odt_2d(pOracle);
    }

	template <class Oracle>
	void one_odt_2d(Oracle* pOracle)
	{
		for (Vertex_iterator
             v = m_mesh->vertices_begin();
             v != m_mesh->vertices_end();
             v++)
		{
			// skip fixed vertices
			if(m_mesh->is_vertex_grounded(v))
				continue;

			Point_2 cc2d = m_mesh->odt_cc_2d(v);
			//Point_3 cc3d = m_mesh->odt_cc_3d(v);
			//Point_2 cc2d = to_2d(cc3d);
            
			const FT z = pOracle->evaluate_height(cc2d);
			Point_3 lifted_centroid(cc2d.x(), cc2d.y(), z);
            
			Point_2 new_point2d = to_2d(lifted_centroid);
            
			if (m_mesh->valid_embedding(v, new_point2d))
				v->point() = lifted_centroid;
		}    
	}

    
	template <class Facet_tree, class Feature_tree>
	void lloyd(Facet_tree& facet_tree, Feature_tree& feature_tree)
	{
		for (Vertex_iterator
             v = m_mesh->vertices_begin();
             v != m_mesh->vertices_end();
             v++)
		{
			bool boundary = m_mesh->is_border_vertex(v);
			Point_2 centroid2d = m_mesh->voronoi_centroid(v, boundary);
			Point_3 centroid3d = to_3d(centroid2d);
            
			// FIXME: centroid2d may be located outside domain.
            
			Point_3 new_point3d;
			if (boundary)
				new_point3d = feature_tree.closest_point(centroid3d);
			else
				new_point3d = facet_tree.closest_point(centroid3d);
			Point_2 new_point2d = to_2d(new_point3d);
            
			if (m_mesh->valid_embedding(v, new_point2d))
			{
				v->set_2d_point(new_point2d);
				stack_delaunay_flip(v);
			}
		}
	}
    
    // auxiliary //
    
    // measure distance between 3D points, using 2D coordinates only
	FT distance_2d(const Point_3& p, const Point_3& q)
	{
		return std::sqrt(CGAL::squared_distance(to_2d(p), to_2d(q)));
	}
    
    // measure distance between 3D points
	FT distance_3d(const Point_3& p, const Point_3& q)
	{
		return std::sqrt(CGAL::squared_distance(p, q));
	}
    
	Point_3 to_3d(const Point_2& p) const
	{
		return Point_3(p.x(), p.y(), 0.0);
	}
    
	Point_2 to_2d(const Point_3& p) const
	{
		return Point_2(p.x(), p.y());
	}
};

#endif
