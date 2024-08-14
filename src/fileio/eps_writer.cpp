#include "fileio/eps_writer.hpp"

#include "algorithm/triangulating.hpp"
#include "datastructure/bbox.hpp"
#include "datastructure/bezier_curve.hpp"
#include "datastructure/guarded_bezier_curve.hpp"
#include "datastructure/triangulation_constraints.hpp"
#include <fstream>
#include <iostream>

namespace bzmsh
{
using std::max;
using std::ofstream;
using std::ostream;
using std::string;

template <typename T>
bool EPSWriter<T>::write(const string& filename, const BezierMesh<T>& mesh, bool resetScaling)
{
    ofstream out(filename);
    if (!out.good() || mesh.allEdges.empty())
    {
        return false;
    }

    BBox<double> bbox(convertVec<T,double>(mesh.allVertices.front()));
    for (const Vec2<T>& vertexT : mesh.allVertices)
    {
        Vec2<double> vertex = convertVec<T,double>(vertexT);
        bbox.widen(vertex);
    }
    for (const Vec2<T>& ctrlptT : mesh.allCtrlpts)
    {
        Vec2<double> ctrlpt = convertVec<T,double>(ctrlptT);
        bbox.widen(ctrlpt);
    }

    double meanDistSq = 0;
    for (const Edge& e : mesh.allEdges)
    {
        double dist = (mesh.allVertices[e.from] - mesh.allVertices[e.to]).length();
        meanDistSq += dist * dist;
    }
    meanDistSq /= mesh.allEdges.size();
    meanDistSq = sqrt(meanDistSq);
    if (resetScaling)
        scalingFactor = 1.0;
    scalingFactor *= 50 / (meanDistSq);

    shiftX = -(bbox.maxX + bbox.minX) / 2;
    shiftY = -(bbox.maxY + bbox.minY) / 2;

    double paddingX = 0.01 * (bbox.maxX - bbox.minX);
    double paddingY = 0.01 * (bbox.maxY - bbox.minY);

    bbox.minX = bbox.minX - paddingX;
    bbox.minY = bbox.minY - paddingY;
    bbox.maxX = bbox.maxX + paddingX;
    bbox.maxY = bbox.maxY + paddingY;

    out << std::setprecision(10);

    writeHeader(out, bbox);
    writeEdges(out, mesh.allEdges, convertVecVec<T, double>(mesh.allVertices), convertVecVec<T, double>(mesh.allCtrlpts), convertStdVec<T, double>(mesh.allCtrlptsWeights));
    if (mesh.ctrlpointFlipped.size() != mesh.allCtrlpts.size())
    {
        mesh.ctrlpointFlipped.resize(mesh.allCtrlpts.size(), false);
    }
    writeCtrlpts(out, convertVecVec<T,double>(mesh.allCtrlpts), mesh.ctrlpointFlipped);
    writeVertices(out, convertVecVec<T,double>(mesh.allVertices));

    out << endl;

    return true;
}

template bool EPSWriter<double>::write(const string& filename, const BezierMesh<double>& mesh, bool resetScaling);
template bool EPSWriter<mpq_class>::write(const string& filename, const BezierMesh<mpq_class>& mesh, bool resetScaling);

template <typename T>
void EPSWriter<T>::writeHeader(ostream& out, const BBox<double>& bbox) const
{
    /// Header
    out << "%!PS-Adobe-3.0 EPSF-3.0"
        << "\n"
        << "%%Pages: 1"
        << "\n"
        << "%%BoundingBox: " << scalingFactor * (bbox.minX + shiftX) << " "
        << scalingFactor * (bbox.minY + shiftY) << " " << scalingFactor * (bbox.maxX + shiftX)
        << " " << scalingFactor * (bbox.maxY + shiftY) << "\n";
    out << "%%EndComments"
        << "\n";

    out << "%%Page: 1 1"
        << "\n";
}

template <typename T>
void EPSWriter<T>::writeEdges(ostream& out,
                           const vector<Edge>& allEdges,
                           const vector<Vec2<double>>& allVertices,
                           const vector<Vec2<double>>& allCtrlpts,
                           const vector<double>& allCtrlptsWeights)
{
    out << "\n"
        << "%Edges section"
        << "\n";
    writeLineStyle(out, colorEdges, lineStyleEdges, lineWidthEdges);

    for (const Edge& e : allEdges)
    {
        if (e.marker == Edge::MARK_DEFAULT || e.marker == Edge::MARK_BOUNDARY
            || e.marker == Edge::MARK_BBOX)
        {
            assert(e.from != -1 && e.to != -1);
            //writeLine(out, LineSegment<double>(allVertices[e.from], allVertices[e.to]));
            vector<Vec2<double>> ctrlpts = {allVertices[e.from]};
            vector<double>weights = { 1.0 };
            for (int ctrlptIndex : e.ctrlpts)
            {
                ctrlpts.emplace_back(allCtrlpts[ctrlptIndex]);
                weights.emplace_back(allCtrlptsWeights[ctrlptIndex]);
            }
            ctrlpts.emplace_back(allVertices[e.to]);
            weights.emplace_back(1.0);
            // Vector graphics only support cubic beziers, other degrees have to be sampled
            writeCurve(out, ctrlpts, weights);
        }
    }

    out << "\n"
        << "%Curves section"
        << "\n";
    writeLineStyle(out, colorCurves, lineStyleCurves, lineWidthCurves);

    for (const Edge& e : allEdges)
    {
        if (e.marker < Edge::MARK_CURVEIDLIMIT)
            continue;
        assert(e.from != -1 && e.to != -1);
        vector<Vec2<double>> ctrlpts = {allVertices[e.from]};
        vector<double>weights = { 1.0 };
        for (int ctrlptIndex : e.ctrlpts)
        {
            ctrlpts.emplace_back(allCtrlpts[ctrlptIndex]);
            weights.emplace_back(allCtrlptsWeights[ctrlptIndex]);
        }
        ctrlpts.emplace_back(allVertices[e.to]);
        weights.emplace_back(1);
        bool isRational = (weights != vector<double>(weights.size(), 1.0));
        if (isRational)
        {
            writeLineStyle(out, colorConic, lineStyleCurves, lineWidthCurves);
            writeCurve(out, ctrlpts, weights);
            writeLineStyle(out, colorCurves, lineStyleCurves, lineWidthCurves);
        }
        else
        {
            writeCurve(out, ctrlpts, weights);
        }
    }
}

template <typename T>
void EPSWriter<T>::writeVertices(ostream& out, const vector<Vec2<double>>& allVertices) const
{
    out << "\n"
        << "%Vertices section"
        << "\n";
    writeLineStyle(out, colorVertices, LineStyle::NONE, 0.0);
    for (const Vec2<double>& vertex : allVertices)
    {
        writeCircle(out, vertex, radiusVertices);
    }
}

template <typename T>
void EPSWriter<T>::writeCtrlpts(ostream& out, const vector<Vec2<double>>& allCtrlpts, const vector<bool>& flipped) const
{
    out << "\n"
        << "%Ctrlpoints section"
        << "\n";
    writeLineStyle(out, colorCtrlpts, LineStyle::NONE, 0.0);
    int i = 0;
    for (const Vec2<double>& ctrlpt : allCtrlpts)
    {
      if(flipped[i])
      {
        writeCircle(out, ctrlpt, radiusCtrlpts*2);
      }
      else
        writeCircle(out, ctrlpt, radiusCtrlpts);
        i++;
    }
}

template <typename T>
void EPSWriter<T>::writeCircle(ostream& out, const Vec2<double>& center, double radius) const
{
    out << (center.x + shiftX) * scalingFactor << " " << (center.y + shiftY) * scalingFactor << " "
        << radius << " 0.0 360.0 arc closepath"
        << "\n"
        << "fill"
        << "\n";
}

template <typename T>
void EPSWriter<T>::writeLineStyle(ostream& out, Color c, LineStyle lineStyle, double lineWidth) const
{
    out << c.r << " " << c.g << " " << c.b << " setrgbcolor"
        << "\n"
        << lineWidth << " setlinewidth"
        << "\n";
    switch (lineStyle)
    {
    case LineStyle::SOLID:
        out << "[" << lineWidth << " 0.0] 0 setdash"
            << "\n";
        break;
    case LineStyle::DASHED:
        out << "[" << 9.0 * lineWidth << " " << 9.0 * lineWidth << "] 0 setdash"
            << "\n";
        break;
    case LineStyle::DOTTED:
        out << "[" << 1.0 * lineWidth << " " << 4.0 * lineWidth << "] 0 setdash"
            << "\n";
        break;
    case LineStyle::DASHDOTTED:
        out << "[" << 9.0 * lineWidth << " " << 4.0 * lineWidth << " " << 1.0 * lineWidth << " "
            << 4.0 * lineWidth << "] 0 setdash"
            << "\n";
        break;
    case LineStyle::NONE:
        out << "[0.0 1.0] 0 setdash"
            << "\n";
        break;
    }
}

template <typename T>
void EPSWriter<T>::writeCurve(ostream& out, const vector<Vec2<double>>& ctrlpts, const vector<double>& weights) const
{
    if (ctrlpts.size() == 2)
    {
        writeLine(out, ctrlpts);
    }
    else if (ctrlpts.size() == 3 && weights == vector<double>(3, 1.0))
    {
        BezierCurve<double> curve(ctrlpts);
        curve.increaseDegree(1);
        writeCubicBezier(out, curve.m_ctrlpts);
    }
    else if (ctrlpts.size() == 3 && weights != vector<double>(3, 1.0))
    {
        writeSampledBezier(out, BezierCurve<double>(ctrlpts, weights), samplesPerCurve);
    }
    else if (ctrlpts.size() == 4 && weights == vector<double>(4, 1.0))
    {
        writeCubicBezier(out, ctrlpts);
    }
    else if (ctrlpts.size() == 4 && weights != vector<double>(4, 1.0))
    {
        writeSampledBezier(out, BezierCurve<double>(ctrlpts, weights), samplesPerCurve);
    }
    else if (ctrlpts.size() >= 5)
    {
        writeSampledBezier(out, BezierCurve<double>(ctrlpts, weights), samplesPerCurve);
    }
}

template <typename T>
void EPSWriter<T>::writeCubicBezier(ostream& out, const vector<Vec2<double>>& c) const
{
    assert(c.size() == 4);
    out << "newpath"
        << "\n"
        << scalingFactor * (c[0].x + shiftX) << " " << scalingFactor * (c[0].y + shiftY)
        << " moveto"
        << "\n"
        << scalingFactor * (c[1].x + shiftX) << " " << scalingFactor * (c[1].y + shiftY) << " "
        << scalingFactor * (c[2].x + shiftX) << " " << scalingFactor * (c[2].y + shiftY) << " "
        << scalingFactor * (c[3].x + shiftX) << " " << scalingFactor * (c[3].y + shiftY)
        << " curveto"
        << "\n"
        << "stroke"
        << "\n";
}

template <typename T>
void EPSWriter<T>::writeSampledBezier(ostream& out,
                                   const BezierCurve<double>& c,
                                   int sampleCount) const
{
    assert(c.degree() >= 0);
    out << "newpath"
        << "\n"
        << scalingFactor * (c[0].x + shiftX) << " " << scalingFactor * (c[0].y + shiftY)
        << " moveto"
        << "\n";
    for (int sample = 1; sample < sampleCount; sample++)
    {
        double t = double(sample) / double(sampleCount);
        Vec2<double> pt = (c(t) + Vec2<double>(shiftX, shiftY)) * scalingFactor;
        out << pt.x << " " << pt.y << " lineto"
            << "\n";
    }
    out << scalingFactor * (c[c.degree()].x + shiftX) << " "
        << scalingFactor * (c[c.degree()].y + shiftY) << " lineto"
        << "\n"
        << "stroke"
        << "\n";
}

template <typename T>
void EPSWriter<T>::writeLine(ostream& out, const LineSegment<double>& l) const
{
    assert(l.degree() >= 0);
    out << "newpath"
        << "\n"
        << scalingFactor * (l.startPt().x + shiftX) << " "
        << scalingFactor * (l.startPt().y + shiftY) << " moveto"
        << "\n"
        << scalingFactor * (l.endPt().x + shiftX) << " " << scalingFactor * (l.endPt().y + shiftY)
        << " lineto"
        << "\n"
        << "stroke"
        << "\n";
}

} // namespace bzmsh
