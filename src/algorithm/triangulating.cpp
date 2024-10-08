#include "algorithm/triangulating.hpp"

#include "algorithm/intersection_tests.hpp"
#include "datastructure/guarded_bezier_curve.hpp"
#include "datastructure/triangulation_constraints.hpp"
#include "fileio/eps_writer.hpp"
#include "util/helper_functions.hpp"
#include "util/logger.hpp"
#include <map>
#include <triangle/triangle.h>

#ifdef ALL_EXACT
#include "datastructure/cgal_types.hpp"
#endif

namespace bzmsh
{

#ifdef ALL_EXACT
bool performCDT(TriangulationConstraints<mpq_class>& in, BezierMesh<mpq_class>& out, bool addBbox)
{
    CDT_2D_EXACT cdt_mesh;

    if (addBbox && !in.addBbox())
        return false;

    std::map<Point_2_exact, int> point_map;
    for( unsigned int i=0; i<in.m_seqPoints.size(); i+=2 )
    {
        //Point_2_exact p( FT_exact( in.m_seqPoints[i] ), FT_exact(in.m_seqPoints[i+1]));
        Point_2_exact p(FT_exact(convertRation(in.m_seqPoints[i])), FT_exact(convertRation(in.m_seqPoints[i + 1])));
        point_map[p] = in.get_point_index(Vec2<mpq_class>(in.m_seqPoints[i], in.m_seqPoints[i+1]));
    }

    std::map< std::pair<Point_2_exact, Point_2_exact>, int> seg_map;
    for( unsigned int i=0; i<in.m_seqSegments.size(); i+=2 )
    {
        int p1_id = 2*in.m_seqSegments[i];
        int p2_id = 2*in.m_seqSegments[i+1];
        //Point_2_exact p1( FT_exact( in.m_seqPoints[p1_id] ), FT_exact(in.m_seqPoints[p1_id+1]));
        //Point_2_exact p2( FT_exact( in.m_seqPoints[p2_id] ), FT_exact(in.m_seqPoints[p2_id+1]));
        Point_2_exact p1(FT_exact(convertRation(in.m_seqPoints[p1_id])), FT_exact(convertRation(in.m_seqPoints[p1_id + 1])));
        Point_2_exact p2(FT_exact(convertRation(in.m_seqPoints[p2_id])), FT_exact(convertRation(in.m_seqPoints[p2_id + 1])));

        seg_map[ std::make_pair(p1, p2) ] = i/2;
        seg_map[ std::make_pair(p2, p1) ] = i/2;

        cdt_mesh.insert_constraint(p1, p2);
    }


    if (in.m_seqPoints.size()/2 != cdt_mesh.number_of_vertices())
    {
        Logger::lout(Logger::DEBUG)
            << "Segments (envelopes) were probably intersecting. Linear triangulation result is incompatible." << endl;
        
        return false;
    }

    //Copy the triangulation
    std::vector<CDT_Edge> edgelist;
    vector<int> edgemarkerlist;
    for(auto e=cdt_mesh.finite_edges_begin(); e != cdt_mesh.finite_edges_end(); ++e)
    {
        edgelist.push_back(*e);

        Point_2_exact p1 = e->first->vertex( (e->second+1)%3 )->point();
        Point_2_exact p2 = e->first->vertex( (e->second+2)%3 )->point();

        if( seg_map.end() != seg_map.find( std::make_pair(p1, p2) ))
            edgemarkerlist.push_back( in.m_seqSegmentMarkers[ seg_map[ std::make_pair(p1, p2) ] ] );
        else
            edgemarkerlist.push_back(0);
    }

    std::vector<CDT_Face> facelist;
    std::set<CDT_Face> boundary_facelist;
    for(auto f=cdt_mesh.finite_faces_begin(); f != cdt_mesh.finite_faces_end(); ++f)
    {
        facelist.push_back(f);
        for(int i=0; i<3; i++)
        {
            if( cdt_mesh.is_infinite( f->neighbor(i)) )
                boundary_facelist.insert(f);
        }
    }

    vector<bool> removeTriangle(cdt_mesh.number_of_faces(), false);
    vector<bool> containsCurve(cdt_mesh.number_of_vertices(), false);
    vector<bool> removeEdge(edgelist.size(), false);

    int add_edge = 0;
    for(unsigned int i=0; i<edgelist.size(); i++)
    {
        if (!removeEdge[i])
        {
            int from = point_map[edgelist[i].first->vertex( (edgelist[i].second+1)%3 )->point()];
            int to = point_map[edgelist[i].first->vertex( (edgelist[i].second+2)%3 )->point()];

            Vec2<mpq_class> vfrom = in.get_point(from);
            Vec2<mpq_class> vto = in.get_point(to);
            out.addEdge(
                vfrom, vto, edgemarkerlist[i] - TriangulationConstraints<mpq_class>::MARK_OFFSET);
            add_edge++;
        }
    }
    int face_added = 0;
    for (int i = 0; i < facelist.size(); i++)
    {
        if (!removeTriangle[i])
        {
            std::vector<int> indices;
            vector<Vec2<mpq_class>> vertices;
            for(int j=0; j<3; j++)
                vertices.push_back( in.get_point( point_map[ facelist[i]->vertex(j)->point() ] ));
            out.addTriangle(vertices);
            face_added++;
        }
    }
    return true;
}
#endif

bool performCDT(TriangulationConstraints<double>& in, BezierMesh<double>& out, bool addBbox)
{
    if (addBbox && !in.addBbox())
        return false;

    // Init what needs to be passed into triangulation code
    triangulateio t_in = {};
    t_in.numberofpoints = in.m_seqPoints.size() / 2;
    if (t_in.numberofpoints > 0)
        t_in.pointlist = in.m_seqPoints.data();
    t_in.numberofsegments = in.m_seqSegments.size() / 2;
    if (t_in.numberofsegments > 0)
    {
        t_in.segmentlist = in.m_seqSegments.data();
        t_in.segmentmarkerlist = in.m_seqSegmentMarkers.data();
    }
    t_in.numberofholes = 0;
    t_in.numberofregions = 0;
    t_in.regionlist = t_in.holelist = nullptr;
    t_in.pointmarkerlist = nullptr;

    // init things with nullptrs that we want in output
    triangulateio t_out = {};
    t_out.pointlist = nullptr;
    t_out.segmentlist = t_out.trianglelist = t_out.edgelist = nullptr;
    t_out.pointmarkerlist = t_out.segmentmarkerlist = nullptr;

    // z for indexing starting from 0, not 1
    // p for specifying segment constraints
    // Q for quiet
    // C for checking consistency
#ifdef NDEBUG
    char flags[] = "zpQe";
#else
    char flags[] = "zpCe";
#endif

    triangulate(flags, &t_in, &t_out, nullptr);

    // Check consistency between input and output
    if (t_out.numberofpoints != t_in.numberofpoints
        || t_out.numberofsegments != t_in.numberofsegments)
    {
        Logger::lout(Logger::DEBUG)
            << "Segments (envelopes) were probably intersecting. Linear triangulation result is incompatible." << endl;
        
        return false;
    }

    // Remove all triangles not enclosed within bbox or envelopes
    vector<bool> removeTriangle(t_out.numberoftriangles, false);
    vector<bool> containsCurve(t_out.numberoftriangles, false);
    vector<bool> removeEdge(t_out.numberofedges, false);

    if (!addBbox)
    {
        map<Edge, int, Edge::compare> edge2index;
        map<int, vector<int>> triangle2edge;

        for (int i = 0; i < t_out.numberofedges; i++)
        {
            edge2index[Edge(t_out.edgelist[2 * i],
                            t_out.edgelist[2 * i + 1],
                            t_out.edgemarkerlist[i] - TriangulationConstraints<double>::MARK_OFFSET)]
                = i;
            if (t_out.edgemarkerlist[i] == 1)
            {
                removeEdge[i] = true;
            }
        }

        for (int i = 0; i < t_out.numberoftriangles; i++)
        {
            vector<int> vertices = {t_out.trianglelist[3 * i],
                                    t_out.trianglelist[3 * i + 1],
                                    t_out.trianglelist[3 * i + 2]};
            vector<int> edges;
            for (int j = 0; j < 3; j++)
            {
                int from = vertices[j];
                int to = vertices[(j + 1) % 3];
                triangle2edge[i].emplace_back(edge2index[Edge(from, to)]);
            }
            for (int edgeindex : triangle2edge[i])
            {
                if (t_out.edgemarkerlist[edgeindex]
                    >= Edge::MARK_CURVEIDLIMIT + TriangulationConstraints<double>::MARK_OFFSET)
                {
                    containsCurve[i] = true;
                }
            }
        }

        bool triangleRemoved = true;
        while (triangleRemoved)
        {
            triangleRemoved = false;
            for (int i = 0; i < t_out.numberoftriangles; i++)
            {
                if (removeTriangle[i])
                    continue;
                bool boundary = false;
                for (int edgeindex : triangle2edge[i])
                {
                    if (removeEdge[edgeindex])
                    {
                        boundary = true;
                    }
                    else if (boundary && t_out.edgemarkerlist[edgeindex] < 2)
                    {
                        removeEdge[edgeindex] = true;
                    }
                }
                if (boundary)
                {
                    removeTriangle[i] = true;
                    triangleRemoved = true;
                }
            }
        }
    }

    // Add into the resulting mesh any mesh elements that have NOT been removed
    for (int i = 0; i < t_out.numberofedges; i++)
    {
        if (!removeEdge[i])
        {
            int from = t_out.edgelist[2 * i];
            int to = t_out.edgelist[2 * i + 1];
            Vec2<double> vfrom
                = Vec2<double>(t_out.pointlist[2 * from], t_out.pointlist[2 * from + 1]);
            Vec2<double> vto = Vec2<double>(t_out.pointlist[2 * to], t_out.pointlist[2 * to + 1]);
            out.addEdge(
                vfrom, vto, t_out.edgemarkerlist[i] - TriangulationConstraints<double>::MARK_OFFSET);
        }
    }

    for (int i = 0; i < t_out.numberoftriangles; i++)
    {
        if (!removeTriangle[i])
        {
            vector<int> indices = {t_out.trianglelist[3 * i],
                                   t_out.trianglelist[3 * i + 1],
                                   t_out.trianglelist[3 * i + 2]};
            vector<Vec2<double>> vertices = {
                Vec2<double>(t_out.pointlist[2 * indices[0]], t_out.pointlist[2 * indices[0] + 1]),
                Vec2<double>(t_out.pointlist[2 * indices[1]], t_out.pointlist[2 * indices[1] + 1]),
                Vec2<double>(t_out.pointlist[2 * indices[2]], t_out.pointlist[2 * indices[2] + 1]),
            };
            out.addTriangle(vertices);
        }
    }

    // Cleanup dynamically allocated resources and return
    trifree(t_out.pointlist);
    trifree(t_out.segmentlist);
    trifree(t_out.trianglelist);
    trifree(t_out.edgelist);
    trifree(t_out.pointmarkerlist);
    trifree(t_out.segmentmarkerlist);
    trifree(t_out.edgemarkerlist);

    return true;
}

template <typename T>
void interpolateInnerCtrlpts(const vector<Vec2<T>>& vertices,
                             int degree,
                             vector<Vec2<T>>& ctrlpts)
{
    // Uniformly distribute control points
    if (degree <= 0)
    {
        ctrlpts.emplace_back(vertices[0]);
        return;
    }
    if (degree == 1)
    {
        ctrlpts.insert(ctrlpts.end(), vertices.begin(), vertices.end());
        return;
    }
    if (degree >= 2)
    {
        ctrlpts.insert(ctrlpts.end(), vertices.begin(), vertices.end());
        for (int index = 0; index < 3; index++)
        {
            LineSegment<T> edge(vertices[index], vertices[(index + 1) % 3]);
            for (int i = 1; i < degree; i++)
            {
                T t = (T)i / (T)degree;
                ctrlpts.emplace_back(edge(t));
            }
        }
        if (degree >= 3)
        {
            vector<Vec2<T>> innerPoints;
            T u = (T)degree - 2 / (T)degree;
            T v = (T)1 / (T)degree;
            T w = (T)1 / (T)degree;
            for (int i = 0; i < 3; i++)
            {
                innerPoints.emplace_back(vertices[0] * u + vertices[1] * v + vertices[2] * w);
                std::swap(u, v); // rotate around
                std::swap(u, w);
            }
            interpolateInnerCtrlpts(innerPoints, degree - 3, ctrlpts);
        }
    }
}
// explicit instantiations
template void interpolateInnerCtrlpts<double>(const vector<Vec2<double>>& vertices,
                             int degree,
                             vector<Vec2<double>>& ctrlpts);
template void interpolateInnerCtrlpts<mpq_class>(const vector<Vec2<mpq_class>>& vertices,
                             int degree,
                             vector<Vec2<mpq_class>>& ctrlpts);


template <typename T>
void calculateInnerCtrlpts(const GuardedBezierCurve<T>& gc, Triangle& tri, BezierMesh<T>& mesh)
{
    int degree = gc.m_curve.degree();

    Vec2<T> p1 = mesh.allVertices[tri.vertices[0]];
    Vec2<T> p2 = mesh.allVertices[tri.vertices[1]];
    Vec2<T> p3 = mesh.allVertices[tri.vertices[2]];
    assert(orientation(p1, p2, p3) == ORI_LEFT);

    Vec2<T> edgeP1P3 = p3 - p1;
    Vec2<T> edgeP2P3 = p3 - p2;

    Vec2<T> slope = gc.m_slope;
    Vec2<T> slope1 = gc.m_slope1;

    vector<Vec2<T>> cpts;
    Vec2<T> hullPt;
    for (const Vec2<T>& ctrlpt : gc.m_curve.m_ctrlpts)
    {
        cpts.emplace_back(ctrlpt);
    }
    assert(orientation(p1, p2, gc.m_guardPtLeft) != ORI_COLINEAR);
    assert(orientation(p1, p2, gc.m_guardPtLeft) != ORI_COLINEAR);
    if (orientation(p1, p2, gc.m_guardPtLeft) == ORI_LEFT)
    {
        hullPt = gc.m_hullPtLeft;
    }
    else
    {
        hullPt = gc.m_hullPtRight;
        std::reverse(cpts.begin(), cpts.end());
    }

    //Find r
    Vec2<T> px = p1;
    if (!getLineIntersection(p2, p3, cpts[ cpts.size()-2], cpts[ cpts.size()-2] + slope1, px))
    {
        Logger::lout(Logger::ERROR)
            << "Regular control points could not be computed for " << gc << endl;
    }
    //Ensure it separates p3 from curve control points
    if(orientation(px, px+slope, hullPt) == orientation(px, px+slope, p3))
        px = hullPt;

    //Get first control point
    Vec2<T> pMin0;
    if (!getLineIntersection(p1, p3, px, px + slope, pMin0))
    {
        Logger::lout(Logger::ERROR)
            << "Regular control points could not be computed for " << gc << endl;
    }
    
    T alpha = 0.2;
    pMin0 = pMin0 * (T(1.0) - alpha) + p3 * alpha;

    Vec2<T> p0 = pMin0;
    if (p1 == p0 || p2 == p0 || p3 == p0)
        p0 = (p1 + p3) * 0.5;

    Vec2<T> pn;
    if (!getLineIntersection(p3, p2, p0, p0 + slope, pn))
    {
        Logger::lout(Logger::ERROR)
            << "Regular control points could not be computed for " << gc << endl;
    }

    // Get the first row of inner control points using cone-line intersections
    vector<Vec2<T>> proj;
    proj.push_back(p0);
    Vec2<T> prev = p0;
    for (uint i = 2; i < cpts.size() - 1; i++)
    {
        Vec2<T> pA, pB;
        if (!getLineIntersection(cpts[i], cpts[i] + edgeP2P3, p0, p0 + slope, pA))
            pA = prev;
        if (!getLineIntersection(cpts[i], cpts[i] + edgeP1P3, p0, p0 + slope, pB))
            pB = pn;

        if (orientation(cpts[i], cpts[i] + edgeP2P3, prev)
                == orientation(cpts[i], cpts[i] + edgeP2P3, pn)
            || orientation(cpts[i], cpts[i] + edgeP2P3, pn) == ORI_COLINEAR)
            pA = prev;
        if (orientation(cpts[i], cpts[i] + edgeP1P3, prev)
                == orientation(cpts[i], cpts[i] + edgeP1P3, pn)
            || orientation(cpts[i], cpts[i] + edgeP1P3, prev) == ORI_COLINEAR)
            pB = pn;

        // final point is average of both extremities
        prev = (pA + pB) * T(0.5);
        assert(prev != hullPt);
        proj.push_back(prev);
    }
    assert(pn != hullPt);
    proj.push_back(pn);

    vector<Vec2<T>> ctrlptsRowWise = proj;

    // Get rest of control points as edge-aligned line intersections from
    // neighboring control points of previous row
    int level = 2;
    while (proj.size() > 2)
    {
        // Get new set of points by intersection of edge paralel lines
        vector<Vec2<T>> projTemp;

        Vec2<T> p_new;
        if (!getLineIntersection(proj[proj.size()-2], proj[proj.size()-2] + edgeP1P3, p2, p3, p_new))
        {
            Logger::lout(Logger::ERROR)
            << "Regular control points could not be computed for " << gc << endl;
            assert(false);
        };

        Vec2<T> p;
        if (!getLineIntersection(p_new, p_new + slope, p1, p3, p))
        {
            Logger::lout(Logger::ERROR)
            << "Regular control points could not be computed for " << gc << endl;
            assert(false);
        };
        projTemp.push_back(p);

        for (uint i = 1; i < proj.size()-2; i++)
        {
            if (!getLineIntersection( proj[i], proj[i] + edgeP1P3, p_new, p_new + slope, p ))
            {
                Logger::lout(Logger::ERROR)
                << "Regular control points could not be computed for " << gc << endl;
                assert(false);
            };
            assert(p != hullPt);
            projTemp.push_back(p);
        }
        projTemp.push_back(p_new);

        proj = projTemp;
        // Fill these new points
        ctrlptsRowWise.insert(ctrlptsRowWise.end(), proj.begin(), proj.end());
        level++;
    }

    // prepend dummy missing row, append dummy missing row
    ctrlptsRowWise.insert(ctrlptsRowWise.begin(), degree + 1, Vec2<T>());
    ctrlptsRowWise.emplace_back(Vec2<T>());
    vector<Vec2<T>> ctrlptsRecursive;
    reorderToRecursive(ctrlptsRowWise, ctrlptsRecursive, degree, false);

    // Extract non-curve-edge control points and inner control points from recursively ordered
    // control points (cf helper_functions.hpp for ordering explanaition)
    vector<Vec2<T>> ctrlptsEdgeP2P3(ctrlptsRecursive.begin() + 3 + 1 * (degree - 1),
                                         ctrlptsRecursive.begin() + 3 + 2 * (degree - 1));
    vector<Vec2<T>> ctrlptsEdgeP3P1(ctrlptsRecursive.begin() + 3 + 2 * (degree - 1),
                                         ctrlptsRecursive.begin() + 3 * degree);
    vector<Vec2<T>> ctrlptsInner(ctrlptsRecursive.begin() + 3 * degree,
                                      ctrlptsRecursive.end());

    vector<int> ctrlptsInnerIndices, ctrlptsEdgeP2P3Indices, ctrlptsEdgeP3P1Indices;
    for (const Vec2<T>& ctrlptInner : ctrlptsInner)
    {
        ctrlptsInnerIndices.emplace_back(mesh.addCtrlpt(ctrlptInner));
    }
    for (const Vec2<T>& ctrlptEdgeP2P3 : ctrlptsEdgeP2P3)
    {
        ctrlptsEdgeP2P3Indices.emplace_back(mesh.addCtrlpt(ctrlptEdgeP2P3));
    }
    for (const Vec2<T>& ctrlptEdgeP3P1 : ctrlptsEdgeP3P1)
    {
        ctrlptsEdgeP3P1Indices.emplace_back(mesh.addCtrlpt(ctrlptEdgeP3P1));
    }
    // if edge in mesh is reversed relative to this triangles order, save ctrlpts in reverse order
    if (tri.vertices[1] == mesh.allEdges[tri.edges[1]].to
        && tri.vertices[2] == mesh.allEdges[tri.edges[1]].from)
    {
        std::reverse(ctrlptsEdgeP2P3Indices.begin(), ctrlptsEdgeP2P3Indices.end());
    }
    if (tri.vertices[2] == mesh.allEdges[tri.edges[2]].to
        && tri.vertices[0] == mesh.allEdges[tri.edges[2]].from)
    {
        std::reverse(ctrlptsEdgeP3P1Indices.begin(), ctrlptsEdgeP3P1Indices.end());
    }
    tri.ctrlpts = ctrlptsInnerIndices;
    mesh.allEdges[tri.edges[1]].ctrlpts = ctrlptsEdgeP2P3Indices;
    mesh.allEdges[tri.edges[2]].ctrlpts = ctrlptsEdgeP3P1Indices;
}

// explicit instantiations
template void calculateInnerCtrlpts<double>(const GuardedBezierCurve<double>& gc, Triangle& tri, BezierMesh<double>& mesh);
template void calculateInnerCtrlpts<mpq_class>(const GuardedBezierCurve<mpq_class>& gc, Triangle& tri, BezierMesh<mpq_class>& mesh);

template <typename T>
void calculateConicInnerCtrlpts(const GuardedBezierCurve<T>& gc, Triangle& tri, BezierMesh<T>& mesh)
{
    int curve_degree = gc.m_curve.degree();
    assert(curve_degree == 2);
    Vec2<T> v00 = mesh.allVertices[tri.vertices[0]];
    Vec2<T> v20 = mesh.allVertices[tri.vertices[1]];
    Vec2<T> v02 = mesh.allVertices[tri.vertices[2]];
    assert(orientation(v00, v20, v02) == ORI_LEFT);
    Vec2<T> slope1 = gc.m_slope1;
    Vec2<T> slope2 = gc.m_slope2;
    Vec2<T> v10 = gc.m_curve.m_ctrlpts[1];
    double w10 = gc.m_curve.m_weights[1];
    Vec2<T>v01, v11;
    if (!getLineIntersection(v10, v10 + slope1, v20, v02, v11))
    {
        v11 = (v20 + v02) / 2.0;
    }
    if (!getLineIntersection(v10, v10 + slope2, v00, v02, v01))
    {
        v01 = (v00 + v02) / 2.0;
    }
    double alpha = 0.2;
    v01 = v01 * (1 - alpha) + v02 * alpha;
    v11 = v11 * (1 - alpha) + v02 * alpha;

    vector<Vec2<T>>ctrlptsRowWise = { v00,v10,v20,v01,v11,v02 };
    vector<T>weightsRowWise = { T(1.0),w10,T(1.0),T(1.0),T(1.0),T(1.0) };
    vector<Vec2<T>>newctrlptsRowWise;
    vector<T>newweightsRowWise;
    IncreaseTriangleDegree(ctrlptsRowWise, weightsRowWise, newctrlptsRowWise, newweightsRowWise, mesh.degree - 2);
    vector<Vec2<T>> ctrlptsRecursive;
    vector<T>weightsRecursive;
    reorderToRecursive(newctrlptsRowWise, ctrlptsRecursive, mesh.degree, false);
    reorderToRecursive(newweightsRowWise, weightsRecursive, mesh.degree, false);
    vector<Vec2<T>> ctrlptsEdgeP2P3(ctrlptsRecursive.begin() + 3 + 1 * (mesh.degree - 1),
        ctrlptsRecursive.begin() + 3 + 2 * (mesh.degree - 1));
    vector<T> weightsEdgeP2P3(weightsRecursive.begin() + 3 + 1 * (mesh.degree - 1),
        weightsRecursive.begin() + 3 + 2 * (mesh.degree - 1));
    vector<Vec2<T>> ctrlptsEdgeP3P1(ctrlptsRecursive.begin() + 3 + 2 * (mesh.degree - 1),
        ctrlptsRecursive.begin() + 3 * mesh.degree);
    vector<T> weightsEdgeP3P1(weightsRecursive.begin() + 3 + 2 * (mesh.degree - 1),
        weightsRecursive.begin() + 3 * mesh.degree);
    vector<Vec2<T>> ctrlptsInner(ctrlptsRecursive.begin() + 3 * mesh.degree,
        ctrlptsRecursive.end());
    vector<T> weightsInner(weightsRecursive.begin() + 3 * mesh.degree,
        weightsRecursive.end());

    vector<int> ctrlptsInnerIndices, ctrlptsEdgeP2P3Indices, ctrlptsEdgeP3P1Indices;
    int k = 0;
    for (const Vec2<T>& ctrlptInner : ctrlptsInner)
    {
        int id = mesh.addCtrlpt(ctrlptInner);
        ctrlptsInnerIndices.emplace_back(id);
        mesh.allCtrlptsWeights[id] = weightsInner[k++];
    }
    k = 0;
    for (const Vec2<T>& ctrlptEdgeP2P3 : ctrlptsEdgeP2P3)
    {
        int id = mesh.addCtrlpt(ctrlptEdgeP2P3);
        ctrlptsEdgeP2P3Indices.emplace_back(id);
        mesh.allCtrlptsWeights[id] = weightsEdgeP2P3[k++];
    }
    k = 0;
    for (const Vec2<T>& ctrlptEdgeP3P1 : ctrlptsEdgeP3P1)
    {
        int id = mesh.addCtrlpt(ctrlptEdgeP3P1);
        ctrlptsEdgeP3P1Indices.emplace_back(id);
        mesh.allCtrlptsWeights[id] = weightsEdgeP3P1[k++];
    }
    if (tri.vertices[1] == mesh.allEdges[tri.edges[1]].to
        && tri.vertices[2] == mesh.allEdges[tri.edges[1]].from)
    {
        std::reverse(ctrlptsEdgeP2P3Indices.begin(), ctrlptsEdgeP2P3Indices.end());
    }
    if (tri.vertices[2] == mesh.allEdges[tri.edges[2]].to
        && tri.vertices[0] == mesh.allEdges[tri.edges[2]].from)
    {
        std::reverse(ctrlptsEdgeP3P1Indices.begin(), ctrlptsEdgeP3P1Indices.end());
    }
    tri.ctrlpts = ctrlptsInnerIndices;
    mesh.allEdges[tri.edges[1]].ctrlpts = ctrlptsEdgeP2P3Indices;
    mesh.allEdges[tri.edges[2]].ctrlpts = ctrlptsEdgeP3P1Indices;
}

// explicit instantiations
template void calculateConicInnerCtrlpts<double>(const GuardedBezierCurve<double>& gc, Triangle& tri, BezierMesh<double>& mesh);
template void calculateConicInnerCtrlpts<mpq_class>(const GuardedBezierCurve<mpq_class>& gc, Triangle& tri, BezierMesh<mpq_class>& mesh);

template <typename T>
void calculatRcubicInnerCtrlpts(const GuardedBezierCurve<T>& gc, Triangle& tri, BezierMesh<T>& mesh)
{
    int curve_degree = gc.m_curve.degree();
    assert(curve_degree == 3);
    Vec2<T> v00, v10, v20, v30, v03;
    v00 = mesh.allVertices[tri.vertices[0]];
    v30 = mesh.allVertices[tri.vertices[1]];
    v03 = mesh.allVertices[tri.vertices[2]];
    assert(orientation(v00, v30, v03) == ORI_LEFT);
    T w10, w20;
    if (orientation(v00, v30, gc.m_guardPtLeft) == ORI_LEFT)
    {
        v10 = gc.m_curve.m_ctrlpts[1];
        v20 = gc.m_curve.m_ctrlpts[2];
        w10 = gc.m_curve.m_weights[1];
        w20 = gc.m_curve.m_weights[2];
    }
    else
    {
        v10 = gc.m_curve.m_ctrlpts[2];
        v20 = gc.m_curve.m_ctrlpts[1];
        w10 = gc.m_curve.m_weights[2];
        w20 = gc.m_curve.m_weights[1];
    }
    Vec2<T> slope1 = gc.m_slope1;
    Vec2<T> slope2 = gc.m_slope2;

    Vec2<T> v01, v02, v21, v12, v11;
    /*Vec2<T> pl1, pl2, pr1, pr2;
    if (!getLineIntersection(v10, v10 + slope2, v00, v03, pl1))
        pl1 = (v00 + v03) * 0.5;
    if (!getLineIntersection(v10, v10 + slope1, v30, v03, pr1))
        pr1 = (v30 + v03) * 0.5;
    if (!getLineIntersection(v20, v20 + slope2, v00, v03, pl2))
        pl2 = (v00 + v03) * 0.5;
    if (!getLineIntersection(v20, v20 + slope1, v30, v03, pr2))
        pr2 = (v30 + v03) * 0.5;
    double alpha = 0.2;
    v01 = (pl1 - v03).length() < (pl2 - v03).length() ? pl1 * (1.0 - alpha) + v03 * alpha : pl2 * (1.0 - alpha) + v03 * alpha;
    v21 = (pr1 - v03).length() < (pr2 - v03).length() ? pr1 * (1.0 - alpha) + v03 * alpha : pr2 * (1.0 - alpha) + v03 * alpha;*/

    Vec2<T>i1 = v00, i2 = v30;
    vector<Vec2<T>>cps{ v00,v10,v20,v30 };
    for (int i = 0; i < 3; ++i)
    {
        Vec2<T>tmp1, tmp2;
        if (!getLineIntersection(cps[i], cps[i + 1], v00, v03, tmp1))
            tmp1 = (v00 + v03) * 0.5;
        if (!getLineIntersection(cps[i], cps[i + 1], v30, v03, tmp2))
            tmp2 = (v30 + v03) * 0.5;
        if ((tmp1 - v03).length() < (i1 - v03).length())
            i1 = tmp1;
        if ((tmp2 - v03).length() < (i2 - v03).length())
            i2 = tmp2;
    }
    double alpha = 0.2;
    v01 = i1 * (1.0 - alpha) + v03 * alpha;
    v21 = i2 * (1.0 - alpha) + v03 * alpha;

    Vec2<T> q1, q2;
    if (!getLineIntersection(v10, v03, v01, v21, q1))
        q1 = v01;
    if (!getLineIntersection(v20, v03, v01, v21, q2))
        q2 = v21;
    v11 = (q1 + q2) * 0.5;

    if (!getLineIntersection(v10, v11, v30, v03, v12))
        v12 = (v21 + v03) * 0.5;
    if (!getLineIntersection(v20, v11, v00, v03, v02))
        v02 = (v01 + v03) * 0.5;

    vector<Vec2<T>>ctrlptsRowWise = { v00,v10,v20,v30,v01,v11,v21,v02,v12,v03 };
    vector<T>weightsRowWise = { T(1.0),w10,w20,T(1.0),T(1.0),(w10 + w20 + T(7.0)) / T(9.0),T(1.0),T(1.0),T(1.0),T(1.0) };
    vector<Vec2<T>>newctrlptsRowWise;
    vector<T>newweightsRowWise;
    IncreaseTriangleDegree(ctrlptsRowWise, weightsRowWise, newctrlptsRowWise, newweightsRowWise, mesh.degree - 3);
    vector<Vec2<T>> ctrlptsRecursive;
    vector<T>weightsRecursive;
    reorderToRecursive(newctrlptsRowWise, ctrlptsRecursive, mesh.degree, false);
    reorderToRecursive(newweightsRowWise, weightsRecursive, mesh.degree, false);
    vector<Vec2<T>> ctrlptsEdgeP1P2(ctrlptsRecursive.begin() + 3,
        ctrlptsRecursive.begin() + 3 + 1 * (mesh.degree - 1));
    vector<T>weightsEdgeP1P2(weightsRecursive.begin() + 3,
        weightsRecursive.begin() + 3 + 1 * (mesh.degree - 1));
    vector<Vec2<T>> ctrlptsEdgeP2P3(ctrlptsRecursive.begin() + 3 + 1 * (mesh.degree - 1),
        ctrlptsRecursive.begin() + 3 + 2 * (mesh.degree - 1));
    vector<T> weightsEdgeP2P3(weightsRecursive.begin() + 3 + 1 * (mesh.degree - 1),
        weightsRecursive.begin() + 3 + 2 * (mesh.degree - 1));
    vector<Vec2<T>> ctrlptsEdgeP3P1(ctrlptsRecursive.begin() + 3 + 2 * (mesh.degree - 1),
        ctrlptsRecursive.begin() + 3 * mesh.degree);
    vector<T> weightsEdgeP3P1(weightsRecursive.begin() + 3 + 2 * (mesh.degree - 1),
        weightsRecursive.begin() + 3 * mesh.degree);
    vector<Vec2<T>> ctrlptsInner(ctrlptsRecursive.begin() + 3 * mesh.degree,
        ctrlptsRecursive.end());
    vector<T> weightsInner(weightsRecursive.begin() + 3 * mesh.degree,
        weightsRecursive.end());

    vector<int> ctrlptsInnerIndices, ctrlptsEdgeP2P3Indices, ctrlptsEdgeP3P1Indices;
    int k = 0;
    for (const Vec2<T>& ctrlptInner : ctrlptsInner)
    {
        int id = mesh.addCtrlpt(ctrlptInner);
        ctrlptsInnerIndices.emplace_back(id);
        mesh.allCtrlptsWeights[id] = weightsInner[k++];
    }
    k = 0;
    for (const Vec2<T>& ctrlptEdgeP2P3 : ctrlptsEdgeP2P3)
    {
        int id = mesh.addCtrlpt(ctrlptEdgeP2P3);
        ctrlptsEdgeP2P3Indices.emplace_back(id);
        mesh.allCtrlptsWeights[id] = weightsEdgeP2P3[k++];
    }
    k = 0;
    for (const Vec2<T>& ctrlptEdgeP3P1 : ctrlptsEdgeP3P1)
    {
        int id = mesh.addCtrlpt(ctrlptEdgeP3P1);
        ctrlptsEdgeP3P1Indices.emplace_back(id);
        mesh.allCtrlptsWeights[id] = weightsEdgeP3P1[k++];
    }
    if (tri.vertices[1] == mesh.allEdges[tri.edges[1]].to
        && tri.vertices[2] == mesh.allEdges[tri.edges[1]].from)
    {
        std::reverse(ctrlptsEdgeP2P3Indices.begin(), ctrlptsEdgeP2P3Indices.end());
    }
    if (tri.vertices[2] == mesh.allEdges[tri.edges[2]].to
        && tri.vertices[0] == mesh.allEdges[tri.edges[2]].from)
    {
        std::reverse(ctrlptsEdgeP3P1Indices.begin(), ctrlptsEdgeP3P1Indices.end());
    }
    tri.ctrlpts = ctrlptsInnerIndices;
    mesh.allEdges[tri.edges[1]].ctrlpts = ctrlptsEdgeP2P3Indices;
    mesh.allEdges[tri.edges[2]].ctrlpts = ctrlptsEdgeP3P1Indices;
}

// explicit instantiations
template void calculatRcubicInnerCtrlpts<double>(const GuardedBezierCurve<double>& gc, Triangle& tri, BezierMesh<double>& mesh);
template void calculatRcubicInnerCtrlpts<mpq_class>(const GuardedBezierCurve<mpq_class>& gc, Triangle& tri, BezierMesh<mpq_class>& mesh);


template <typename T>
void calculateStraightTCtrlpts(Triangle& tri, BezierMesh<T>& mesh)
{
    int degree = mesh.degree;
    assert(!mesh.allEdges[tri.edges[1]].ctrlpts.empty());
    assert(!mesh.allEdges[tri.edges[2]].ctrlpts.empty());
    assert(mesh.allEdges[tri.edges[0]].ctrlpts.empty());
    Vec2<T> p1 = mesh.allVertices[tri.vertices[0]];
    Vec2<T> p2 = mesh.allVertices[tri.vertices[1]];
    Vec2<T> p3 = mesh.allVertices[tri.vertices[2]];
    vector<Vec2<T>>cpts1;
    for (auto p : mesh.allEdges[tri.edges[1]].ctrlpts)
    {
        cpts1.emplace_back(mesh.allCtrlpts[p]);
    }
    if (tri.vertices[1] != mesh.allEdges[tri.edges[1]].from)
    {
        std::reverse(cpts1.begin(), cpts1.end());
    }
    assert(cpts1.size() == (degree - 1));
    vector<Vec2<T>>cps;
    for (auto p : cpts1)
    {
        Vec2<T> intersect;
        getLineIntersection(p1, p2, p, p + p1 - p3, intersect);
        cps.emplace_back(intersect);
    }
    std::reverse(cps.begin(), cps.end());

    vector<Vec2<T>>innerpoints;
    vector<LineSegment<T>>lines;
    for (int i = 0; i < cps.size(); ++i)
    {
        lines.emplace_back(cps[i], cpts1[cps.size() - i - 1]);
    }
    for (int i = 0; i < lines.size() - 1; ++i)
    {
        for (int j = 0; j < lines.size() - 1 - i; ++j)
        {
            T t = (i + 1) / (T)(degree - 2 - j + 1);
            innerpoints.emplace_back(lines[j](t));
        }
    }
    for (auto p : innerpoints)
    {
        tri.ctrlpts.emplace_back(mesh.addCtrlpt(p));
    }

    if (tri.vertices[0] == mesh.allEdges[tri.edges[0]].to)
    {
        std::reverse(cps.begin(), cps.end());
    }
    for (auto p : cps)
    {
        int id = mesh.addCtrlpt(p);
        mesh.allEdges[tri.edges[0]].ctrlpts.emplace_back(id);
    }
}

// explicit instantiations
template void calculateStraightTCtrlpts<double>(Triangle& tri, BezierMesh<double>& mesh);
template void calculateStraightTCtrlpts<mpq_class>(Triangle& tri, BezierMesh<mpq_class>& mesh);

template <typename T>
bool insertMeshCtrlpts(const vector<GuardedBezierCurve<T>>& guardedCurves, BezierMesh<T>& mesh)
{
    if (guardedCurves.empty())
        return false;
    int degree = mesh.degree;

    // Insert curve controlpoints for all edges that represent curves
    // This does not include first and last control points, those are stored as vertices
    for (Edge& e : mesh.allEdges)
    {
        if (e.marker >= Edge::MARK_CURVEIDLIMIT)
        {
            int id = e.marker;
            BezierCurve<T> curve = guardedCurves[e.marker].m_curve;
            assert(id == curve.m_id);
            if (curve.degree() == 1)
                continue;
            if (curve.isRconic() || curve.isRcubic())
                curve.increaseDegree(mesh.degree - curve.degree());
            for (int i = 1; i < curve.degree(); i++)
            {
                int cp_id = mesh.addCtrlpt(curve[i]);
                e.ctrlpts.emplace_back(cp_id);
                mesh.allCtrlptsWeights[cp_id] = curve.m_weights[i];
            }
            if (orientation(mesh.allVertices[e.from],
                mesh.allVertices[e.to],
                guardedCurves[id].m_guardPtLeft)
                == ORI_RIGHT)
            {
                std::reverse(e.ctrlpts.begin(), e.ctrlpts.end());
            }
        }
    }

    // Calculate inner control points for all triangles.
    // For those not containing a curve this will uniformly distribute the points
    // For those containing a curve this will properly calculate the points to assure
    // injectivity
    if (degree >= 2)
    {
        for (Triangle& t : mesh.allTriangles)
        {
            assert(t.ctrlpts.empty());
            assert(t.edges.size() == 3);
            assert(t.vertices.size() == 3);
            int curveEdge = -1;
            for (int i = 0; i < 3; i++)
            {
                if (mesh.allEdges[t.edges[i]].marker >= Edge::MARK_CURVEIDLIMIT)
                {
                    curveEdge = i;
                    break;
                }
            }
            if (curveEdge != -1)
            {
                // Reorder, so curve edge is always FIRST edge of any triangle
                if (curveEdge == 1)
                {
                    std::swap(t.edges[1], t.edges[0]);
                    std::swap(t.vertices[1], t.vertices[0]);
                    std::swap(t.edges[2], t.edges[1]);
                    std::swap(t.vertices[2], t.vertices[1]);
                }
                else if (curveEdge == 2)
                {
                    std::swap(t.edges[2], t.edges[0]);
                    std::swap(t.vertices[2], t.vertices[0]);
                    std::swap(t.edges[1], t.edges[2]);
                    std::swap(t.vertices[1], t.vertices[2]);
                }
                assert(mesh.allEdges[t.edges[0]].marker >= Edge::MARK_CURVEIDLIMIT);
                assert(mesh.allEdges[t.edges[1]].marker < Edge::MARK_CURVEIDLIMIT);
                assert(mesh.allEdges[t.edges[2]].marker < Edge::MARK_CURVEIDLIMIT);
                GuardedBezierCurve<T>curve_ = guardedCurves[mesh.allEdges[t.edges[0]].marker];
                if (curve_.m_curve.degree() == 1)
                    continue;
                if (curve_.m_curve.isRconic())
                {
                    calculateConicInnerCtrlpts(guardedCurves[mesh.allEdges[t.edges[0]].marker], t, mesh);
                }
                else if (curve_.m_curve.isRcubic())
                {
                    calculatRcubicInnerCtrlpts(guardedCurves[mesh.allEdges[t.edges[0]].marker], t, mesh);
                }
                else
                {
                    calculateInnerCtrlpts(guardedCurves[mesh.allEdges[t.edges[0]].marker], t, mesh);
                }
            }
        }
    }

    if (degree >= 3)
    {
        for (Triangle& t : mesh.allTriangles)
        {
            int curveEdge = -1;
            for (int i = 0; i < 3; i++)
            {
                if (mesh.allEdges[t.edges[i]].marker >= Edge::MARK_CURVEIDLIMIT)
                {
                    curveEdge = i;
                    break;
                }
            }
            if (curveEdge != -1)
                continue;
            std::set<int>curvewithpoints;
            int num = 0;
            for (int i = 0; i < 3; i++)
            {
                if (!mesh.allEdges[t.edges[i]].ctrlpts.empty())
                {
                    curvewithpoints.insert(i);
                }
            }
            if (curvewithpoints.size() == 2)
            {
                if (curvewithpoints.find(1) == curvewithpoints.end())
                {
                    std::swap(t.edges[1], t.edges[0]);
                    std::swap(t.vertices[1], t.vertices[0]);
                    std::swap(t.edges[2], t.edges[1]);
                    std::swap(t.vertices[2], t.vertices[1]);
                }
                else if (curvewithpoints.find(2) == curvewithpoints.end())
                {
                    std::swap(t.edges[2], t.edges[0]);
                    std::swap(t.vertices[2], t.vertices[0]);
                    std::swap(t.edges[1], t.edges[2]);
                    std::swap(t.vertices[1], t.vertices[2]);
                }
                calculateStraightTCtrlpts(t, mesh);
            }
        }
    }

    // For any edges that have not received control points so far, uniformly distribute control
    // points over line segment
    for (Edge& e : mesh.allEdges)
    {
        assert(e.from != -1 && e.to != -1);
        if (e.ctrlpts.size() != degree-1)
        {
            e.ctrlpts.clear();
            LineSegment<T> line(mesh.allVertices[e.from], mesh.allVertices[e.to]);
            for (int i = 1; i < degree; i++)
            {
                T t = (T)i / (T)degree;
                e.ctrlpts.emplace_back(mesh.addCtrlpt(line(t)));
            }
        }
    }

    if (degree >= 3)
    {
        for (Triangle& t : mesh.allTriangles)
        {
            if (t.ctrlpts.size() == 0)
            {
                vector<Vec2<T>> innerCtrlpts;
                vector<Vec2<T>> vertices = { mesh.allVertices[t.vertices[0]],
                                                 mesh.allVertices[t.vertices[1]],
                                                 mesh.allVertices[t.vertices[2]] };
                vector<Vec2<T>> innerPoints;
                T u = ((T)degree - 2) / (T)degree;
                T v = (T)1 / (T)degree;
                T w = (T)1 / (T)degree;
                for (int i = 0; i < 3; i++)
                {
                    innerPoints.emplace_back(vertices[0] * u + vertices[1] * v + vertices[2] * w);
                    std::swap(u, v); // rotate around
                    std::swap(u, w);
                }
                interpolateInnerCtrlpts(innerPoints, degree - 3, innerCtrlpts);
                for (const Vec2<T> ctrlpt : innerCtrlpts)
                {
                    t.ctrlpts.emplace_back(mesh.addCtrlpt(ctrlpt));
                }
            }
        }
    }

    mesh.allCtrlptsWeights.resize(mesh.allCtrlpts.size(), T(1.0));
    mesh.allVerticesWeights.resize(mesh.allVertices.size(), T(1.0));
    return true;
}

// explicit instantiations
template bool insertMeshCtrlpts<double>(const vector<GuardedBezierCurve<double>>& guardedCurves, BezierMesh<double>& mesh);
template bool insertMeshCtrlpts<mpq_class>(const vector<GuardedBezierCurve<mpq_class>>& guardedCurves, BezierMesh<mpq_class>& mesh);

template <typename T>
void IncreaseTriangleDegree(const std::vector<Vec2<T>>& rowwisectlps,
    const std::vector<T>& rowwiseweights,
    std::vector<Vec2<T>>& resultctlps,
    std::vector<T>& resultweights,
    int delta)
{
    if (delta <= 0)
    {
        resultctlps = rowwisectlps;
        resultweights = rowwiseweights;
        return;
    }
    resultctlps.clear();
    resultweights.clear();
    std::vector<std::vector<Vec2<T>>>newctlps, prevctlps;
    std::vector<std::vector<T>>newweights, prevweights;
    if (rowwisectlps.size() == 6)
    {
        prevctlps = { {rowwisectlps[0],rowwisectlps[1],rowwisectlps[2]},{rowwisectlps[3],rowwisectlps[4]},{rowwisectlps[5]} };
        prevweights = { {rowwiseweights[0],rowwiseweights[1],rowwiseweights[2]},{rowwiseweights[3],rowwiseweights[4]},{rowwiseweights[5]} };
    }
    else if (rowwisectlps.size() == 10)
    {
        prevctlps = { {rowwisectlps[0],rowwisectlps[1],rowwisectlps[2],rowwisectlps[3]},{rowwisectlps[4],rowwisectlps[5],rowwisectlps[6]},{rowwisectlps[7],rowwisectlps[8]},{rowwisectlps[9]} };
        prevweights = { {rowwiseweights[0],rowwiseweights[1],rowwiseweights[2],rowwiseweights[3]},{rowwiseweights[4],rowwiseweights[5],rowwiseweights[6]},{rowwiseweights[7],rowwiseweights[8]},{rowwiseweights[9]} };
    }
    while (delta-- > 0)
    {
        IncreaseTriangleOneDegree(prevctlps, prevweights, newctlps, newweights);
        prevctlps = newctlps;
        prevweights = newweights;
    }
    for (const auto& e1 : newctlps)
    {
        for (const auto& e2 : e1)
        {
            resultctlps.emplace_back(e2);
        }
    }
    for (const auto& e1 : newweights)
    {
        for (const auto& e2 : e1)
        {
            resultweights.emplace_back(e2);
        }
    }
}

// explicit instantiations
template void IncreaseTriangleDegree<double>(const std::vector<Vec2<double>>& rowwisectlps,
    const std::vector<double>& rowwiseweights,
    std::vector<Vec2<double>>& resultctlps,
    std::vector<double>& resultweights,
    int delta);
template void IncreaseTriangleDegree<mpq_class>(const std::vector<Vec2<mpq_class>>& rowwisectlps,
    const std::vector<mpq_class>& rowwiseweights,
    std::vector<Vec2<mpq_class>>& resultctlps,
    std::vector<mpq_class>& resultweights,
    int delta);

template <typename T>
void IncreaseTriangleOneDegree(const std::vector<std::vector<Vec2<T>>>& prevctlps,
    const std::vector<std::vector<T>>& prevweights,
    std::vector<std::vector<Vec2<T>>>& newctlps,
    std::vector<std::vector<T>>& newweights)
{
    assert(prevctlps.size() == prevweights.size());
    unsigned int newdegree = prevctlps.size();
    newctlps.resize(newdegree + 1);
    newweights.resize(newdegree + 1);
    assert(newctlps.size() == newweights.size());
    for (unsigned int i = 0; i < newdegree; ++i)
    {
        newctlps[i].resize(prevctlps[i].size() + 1);
        newweights[i].resize(prevweights[i].size() + 1);
    }
    newctlps.back().resize(1);
    newweights.back().resize(1);
    for (int i = 0; i < newctlps.size(); ++i)
    {
        for (int j = 0; j < newctlps[i].size(); ++j)
        {
            Vec2<T>b1, b2, b3;
            T w1 = T(0.0), w2 = w1, w3 = w1;
            if (i - 1 >= 0 && j < prevweights[i - 1].size())
            {
                b1 = prevctlps[i - 1][j] * prevweights[i - 1][j];
                w1 = prevweights[i - 1][j];
            }
            if (i < prevctlps.size() && j - 1 >= 0 && j - 1 < prevctlps[i].size())
            {
                b2 = prevctlps[i][j - 1] * prevweights[i][j - 1];
                w2 = prevweights[i][j - 1];
            }
            if (i < prevctlps.size() && j < prevctlps[i].size() && i + j <= newdegree)
            {
                b3 = prevctlps[i][j] * prevweights[i][j];
                w3 = prevweights[i][j];
            }
            newctlps[i][j] = (b1 * T(i) + b2 * T(j) + b3 * T(newdegree - i - j)) / T(newdegree);
            newweights[i][j] = (w1 * T(i) + w2 * T(j) + w3 * T(newdegree - i - j)) / T(newdegree);
            newctlps[i][j] = newctlps[i][j] / newweights[i][j];
        }
    }
}

// explicit instantiations
template void IncreaseTriangleOneDegree<double>(const std::vector<std::vector<Vec2<double>>>& prevctlps,
    const std::vector<std::vector<double>>& prevweights,
    std::vector<std::vector<Vec2<double>>>& newctlps,
    std::vector<std::vector<double>>& newweights);
template void IncreaseTriangleOneDegree<mpq_class>(const std::vector<std::vector<Vec2<mpq_class>>>& prevctlps,
    const std::vector<std::vector<mpq_class>>& prevweights,
    std::vector<std::vector<Vec2<mpq_class>>>& newctlps,
    std::vector<std::vector<mpq_class>>& newweights);

template <typename T>
void subdivide(const BezierMesh<T>& mesh, BezierMesh<T>& subTriangulation)
{
    if (mesh.allEdges.empty())
        return;

    int degree = mesh.degree;

    for (const Triangle& t : mesh.allTriangles)
    {
        vector<Vec2<T>> ctrlpts;
        for (int vindex : t.vertices)
        {
            ctrlpts.emplace_back(mesh.allVertices[vindex]);
        }
        for (int i = 0; i < 3; i++)
        {
            int eindex = t.edges[i];
            int vfrom = t.vertices[i];
            vector<Vec2<T>> edgeCtrlpts;
            for (int ctrlindex : mesh.allEdges[eindex].ctrlpts)
            {
                edgeCtrlpts.emplace_back(mesh.allCtrlpts[ctrlindex]);
            }
            if (vfrom != mesh.allEdges[eindex].from)
            {
                std::reverse(edgeCtrlpts.begin(), edgeCtrlpts.end());
            }
            ctrlpts.insert(ctrlpts.end(), edgeCtrlpts.begin(), edgeCtrlpts.end());
        }
        for (int innerctrlptindex : t.ctrlpts)
        {
            ctrlpts.emplace_back(mesh.allCtrlpts[innerctrlptindex]);
        }
        assert((int)ctrlpts.size() == (degree + 1) * (degree + 2) / 2);

        vector<int> indC;
        int val = 0;
        int add = degree + 1;
        for (int i = 0; i <= degree; i++)
        {
            indC.push_back(val);
            val += add;
            add--;
        }
        vector<Vec2<T>> ctrlptsRowwise;
        reorderToRowwise(ctrlpts, ctrlptsRowwise, degree);

        vector<vector<Vec2<T>>> subTris(4);
        for (int i = 0; i <= degree; i++)
            for (int j = 0; j <= degree - i; j++)
                subTris[0].push_back(interpolate(ctrlptsRowwise, indC, i, j, i, j, degree - i - j));
        for (int i = 0; i <= degree; i++)
            for (int m = 0; m <= i; m++)
                subTris[1].push_back(interpolate(ctrlptsRowwise, indC, i, 0, m, 0, degree - i));
        for (int i = 0; i <= degree; i++)
            for (int m = 0; m <= degree - i; m++)
                subTris[2].push_back(interpolate(ctrlptsRowwise, indC, i, degree - i, i, m, 0));
        for (int i = 0; i <= degree; i++)
            for (int m = 0; m <= degree - i; m++)
                subTris[3].push_back(interpolate(ctrlptsRowwise, indC, 0, i, 0, i, m));
        std::reverse(subTris[1].begin(), subTris[1].end());

        for (int tindex = 0; tindex < 4; tindex++)
        {
            vector<Vec2<T>> ctrlptsrecursive;
            bool reverse
                = orientation(subTris[tindex][0], subTris[tindex][degree], subTris[tindex].back())
                  == ORI_RIGHT;
            reorderToRecursive(subTris[tindex], ctrlptsrecursive, degree, reverse);
            subTris[tindex] = ctrlptsrecursive;
            vector<int> vs(3), es(3);
            for (int i = 0; i < 3; i++)
            {
                vs[i] = subTriangulation.addVertex(subTris[tindex][i]);
            }
            for (int i = 0; i < 3; i++)
            {
                es[i] = subTriangulation.addEdge(vs[i], vs[(i + 1) % 3], tindex == 0 ? -3 : -4);

                Edge& e = subTriangulation.allEdges[es[i]];
                if (e.ctrlpts.empty())
                {
                    for (int j = 0; j < degree - 1; j++)
                    {
                        e.ctrlpts.emplace_back(
                            subTriangulation.addCtrlpt(subTris[tindex][3 + (degree - 1) * i + j]));
                    }
                    if (e.from != vs[i])
                    {
                        assert(e.to == vs[i]);
                        std::reverse(e.ctrlpts.begin(), e.ctrlpts.end());
                    }
                }
            }
            int t = subTriangulation.addTriangle(vs);
            for (uint i = 3 * degree; i < subTris[tindex].size(); i++)
            {
                subTriangulation.allTriangles[t].ctrlpts.emplace_back(
                    subTriangulation.addCtrlpt(subTris[tindex][i]));
            }
        }
    }
}

// explicit instantiations
template void subdivide<double>(const BezierMesh<double>& mesh, BezierMesh<double>& subTriangulation);
template void subdivide<mpq_class>(const BezierMesh<mpq_class>& mesh, BezierMesh<mpq_class>& subTriangulation);


} // namespace bzmsh
