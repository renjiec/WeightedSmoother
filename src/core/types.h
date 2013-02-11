#ifndef _TYPES_
#define _TYPES_

#undef min
#undef max
#define EPS 1.0e-12

//#include <CGAL/Cartesian.h>
//typedef CGAL::Cartesian<double> Kernel;

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef Kernel::FT FT;
typedef Kernel::Ray_2 Ray_2;
typedef Kernel::Ray_3 Ray_3;
typedef Kernel::Line_2 Line_2;
typedef Kernel::Line_3 Line_3;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Vector_2 Vector_2;
typedef Kernel::Vector_3 Vector_3;
typedef Kernel::Segment_2  Segment_2;
typedef Kernel::Segment_3  Segment_3;
typedef Kernel::Triangle_2 Triangle_2;
typedef Kernel::Triangle_3 Triangle_3;
typedef Kernel::Aff_transformation_3 Aff_transformation;

#include "polygon.h"
typedef CPolygon<Kernel> Polytope;

#include <CGAL/Polygon_2.h>
typedef CGAL::Polygon_2<Kernel> Polygon_2;

#include <CGAL/Bbox_2.h>
typedef CGAL::Bbox_2 Bbox_2;

#include "cdt.h"
#include <CGAL/intersections.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
typedef CGAL::Exact_intersections_tag Itag;
typedef CGAL::Triangulation_vertex_base_2<Kernel> CDT_Vb;
typedef CGAL::Constrained_triangulation_face_base_2<Kernel> CDT_Fb;
typedef CDT_vertex_base<Kernel, CDT_Vb> CDT_CVb;
typedef CDT_face_base<Kernel, CDT_Fb>  CDT_CFb;
typedef CGAL::Triangulation_data_structure_2<CDT_CVb, CDT_CFb> CDT_CTDS;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, CDT_CTDS, Itag> MCDT;
typedef CCDT<MCDT> CDT;

#include "halfsphere.h"
#include "height_cdt.h"
typedef CHeight_cdt<CDT> Height_cdt;

#include "oracle.h"
typedef COracle<Polytope, CDT> Oracle;
typedef Oracle::Domain         Domain;
typedef Oracle::LoadField      LoadField;

#include "variable_load.h"
typedef CVariable_load<Kernel> Variable_load;


#include "mesh.h"
typedef CMesh<Kernel, Enriched_items> Mesh;

#include "builder_ifs.h"
typedef CBuilder_ifs<Mesh> Builder_ifs;

#include "predefined.h"
typedef CBuilder_predefined<Mesh> Builder_predefined;

#include "meshing.h"
typedef CMeshing<Mesh> Meshing;

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>

#include <CGAL/AABB_triangle_primitive.h>
typedef std::list<Triangle_3>::iterator TIterator;
typedef CGAL::AABB_triangle_primitive<Kernel, TIterator> Triangle_primitive;
typedef CGAL::AABB_traits<Kernel, Triangle_primitive> Facet_traits;
typedef CGAL::AABB_tree<Facet_traits> Facet_tree;

#include <CGAL/AABB_segment_primitive.h>
typedef std::list<Segment_3>::iterator SIterator;
typedef CGAL::AABB_segment_primitive<Kernel, SIterator> Segment_primitive;
typedef CGAL::AABB_traits<Kernel, Segment_primitive> Edge_traits;
typedef CGAL::AABB_tree<Edge_traits> Feature_tree;

#include "energy.h"
typedef CEnergy<Meshing, Oracle> Energy;

#include "../third/line_search.h"
#include "../third/newton.h"
#include "../third/lbfgs.h"

#ifdef USE_EIGEN
#include "../third/eigen_cg.h"
typedef EigenLinearSolver LinearSolver;
#else
#include "../third/suite_sparse_qr.h"
typedef SuiteSparseQRFactorizer LinearSolver;
#endif

#include "qenergy.h"
typedef CQEnergy<Energy> QEnergy;
typedef LineSearch<QEnergy> QLineSearch;
typedef LBFGSSolver<QLineSearch> QLBFGSSolver;

#include "wenergy.h"
typedef CWEnergy<Energy, LinearSolver, Facet_tree> WEnergy;
typedef LineSearch<WEnergy> WLineSearch;
typedef NewtonSolver<WLineSearch, LinearSolver> WNewtonSolver;

#include "wipopt.h"
typedef CWIpopt<Energy> WIpopt;

#include "hipopt.h"
typedef CHIpopt<Energy> HIpopt;

#include "eipopt.h"
typedef CEIpopt<Energy> EIpopt;

#include "zopt.h"
typedef CZOpt<Mesh, LinearSolver> ZOpt;

#include "harmonic.h"
typedef CHarmonicSolver<Mesh, LinearSolver> HarmonicSolver;

#include "transparent.h"
typedef Render_transparent<Kernel>::Projected_facet Proj_face;

#endif
