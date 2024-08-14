#ifndef BZMSH_PARTITION_CURVE
#define BZMSH_PARTITION_CURVE

#include "datastructure/bezier_curve.hpp"
#include "datastructure/line_segment.hpp"
#include "util/helper_functions.hpp"

namespace bzmsh
{

/**
 * @brief Reparametrizes a given bezier curve with t_start > 0 or t_end < 1
 *        to t_start = 0, t_end = 1. Original t_start and t_end (referring
 *        to an original curve this curve may have been created from)
 *        are left untouched
 *
 * @tparam T float-type (used for double and mpq_class)
 * @param c curve to crop
 * @return false if ill-defined t_start, t_end were encountered
 * @return true else
 */
template <typename T> bool crop(BezierCurve<T>& c)
{
    if (c.m_tEnd > T(1) || c.m_tEnd <= c.m_tStart || c.m_tStart < T(0))
        return false;

    T vtStart = c.m_origTStart;
    T vtEnd = c.m_origTEnd;
    T tStartNew = (c.m_tEnd == T(1)) ? c.m_tStart : c.m_tStart / c.m_tEnd;
    vector<BezierCurve<T>> bisectedC;
    if (bisect(c, c.m_tEnd, bisectedC))
    {
        c = bisectedC[0];
    }
    if (bisect(c, tStartNew, bisectedC))
    {
        c = bisectedC[1];
    }
    c.m_origTStart = vtStart;
    c.m_origTEnd = vtEnd;
    return true;
}

/**
 * @brief Reparametrizes a given line segment with t_start > 0 or t_end < 1
 *        to t_start = 0, t_end = 1. virtual t_start and t_end (referring
 *        to an original curve this curve may have been created from)
 *        are left untouched
 *
 * @tparam T float-type (used for double and mpq_class)
 * @param l line to crop
 * @return false if ill-defined t_start, t_end were encountered
 * @return true else
 */
template <typename T> bool crop(LineSegment<T>& l)
{
    if (l.m_tEnd > T(1) || l.m_tEnd <= l.m_tStart || l.m_tStart < T(0))
        return false;

    vector<Vec2<T>> ctrlpts_new(2);
    ctrlpts_new[0] = l(l.m_tStart);
    ctrlpts_new[1] = l(l.m_tEnd);
    l.m_tStart = 0;
    l.m_tEnd = 1;
    l.m_ctrlpts = ctrlpts_new;
    return true;
}

/**
 * @brief Partition a given bezier curve into two parts at @p t, yielding
 *        two new curves.
 *
 * @tparam T float-type (used for double and mpq_class)
 * @param c curve to bisect
 * @param t parameter value to bisect at (yielding pieces corresponding to [0, t], [t, 1])
 * @param curvesOut resulting two curve pieces (output)
 * @return false if t outside (0, 1), no bisection performed
 * @return true else
 */
template <typename T>
bool bisect(const BezierCurve<T>& c, const T& t, vector<BezierCurve<T>>& curvesOut)
{
    curvesOut.clear();
    if (t >= T(1) || t <= T(0))
    {
        return false;
    }

    int d = c.degree();
    vector<Vec2<T>> ctrlptsFront, ctrlptsBack;
    ctrlptsFront.emplace_back(c.m_ctrlpts[0]);
    ctrlptsBack.emplace_back(c.m_ctrlpts[d]);

    if (!c.isRconic() && !c.isRcubic())
    {
        std::vector<Vec2<T>> ctrlpts = c.m_ctrlpts;
        for (int i = 0; i < d; i++)
        {
            std::vector<Vec2<T>> ctrlptsNew;
            for (int j = 0; j < d - i; j++)
                ctrlptsNew.emplace_back((ctrlpts[j] * (T(1.0) - t) + ctrlpts[j + 1] * t));

            ctrlpts = ctrlptsNew;
            ctrlptsFront.emplace_back(ctrlpts[0]);
            ctrlptsBack.emplace_back(ctrlpts[ctrlpts.size() - 1]);
        }
        std::reverse(ctrlptsBack.begin(), ctrlptsBack.end());
        ctrlptsBack.front() = ctrlptsFront.back();

        curvesOut
            = { BezierCurve<T>(ctrlptsFront, 0, 1, c.m_id), BezierCurve<T>(ctrlptsBack, 0, 1, c.m_id) };
    }
    else if (c.isRconic())
    {
        Vec2<T>p1 = c.m_ctrlpts[0] * c.m_weights[0] * (T(1.0) - t) + c.m_ctrlpts[1] * c.m_weights[1] * t;
        T w1 = c.m_weights[0] * (T(1.0) - t) + c.m_weights[1] * t;
        p1 = p1 / w1;
        Vec2<T>p2 = c.m_ctrlpts[1] * c.m_weights[1] * (T(1.0) - t) + c.m_ctrlpts[2] * c.m_weights[2] * t;
        T w2 = c.m_weights[1] * (T(1.0) - t) + c.m_weights[2] * t;
        p2 = p2 / w2;
        Vec2<T>p3 = p1 *w1* (T(1.0) - t) + p2* w2* t;
        T w3 = w1 * (T(1.0) - t) + w2 * t;
        p3 = p3 / w3;
        ctrlptsFront.emplace_back(p1);
        ctrlptsFront.emplace_back(p3);
        ctrlptsBack.emplace_back(p2);
        ctrlptsBack.emplace_back(p3);
        std::reverse(ctrlptsBack.begin(), ctrlptsBack.end());
        curvesOut
            = { BezierCurve<T>(ctrlptsFront,{c.m_weights[0],w1,w3}, 0, 1, c.m_id), BezierCurve<T>(ctrlptsBack,{w3,w2,c.m_weights[2]}, 0, 1, c.m_id) };
    }
    else
    {
        Vec2<T>p1 = c.m_ctrlpts[0] * c.m_weights[0] * (T(1.0) - t) + c.m_ctrlpts[1] * c.m_weights[1] * t;
        T w1 = c.m_weights[0] * (T(1.0) - t) + c.m_weights[1] * t;
        p1 = p1 / w1;
        Vec2<T>p2 = c.m_ctrlpts[1] * c.m_weights[1] * (T(1.0) - t) + c.m_ctrlpts[2] * c.m_weights[2] * t;
        T w2 = c.m_weights[1] * (T(1.0) - t) + c.m_weights[2] * t;
        p2 = p2 / w2;
        Vec2<T>p3 = c.m_ctrlpts[2] * c.m_weights[2] * (T(1.0) - t) + c.m_ctrlpts[3] * c.m_weights[3] * t;
        T w3 = c.m_weights[2] * (T(1.0) - t) + c.m_weights[3] * t;
        p3 = p3 / w3;
        ctrlptsFront.emplace_back(p1);
        ctrlptsBack.emplace_back(p3);
        Vec2<T>p4 = p1* w1* (T(1.0) - t) + p2* w2* t;
        T w4 = w1 * (T(1.0) - t) + w2 * t;
        p4 = p4 / w4;
        Vec2<T>p5 = p2* w2* (T(1.0) - t) + p3* w3* t;
        T w5 = w2 * (T(1.0) - t) + w3 * t;
        p5 = p5 / w5;
        ctrlptsFront.emplace_back(p4);
        ctrlptsBack.emplace_back(p5);
        Vec2<T>p6 = p4 *w4* (T(1.0) - t) + p5 *w5* t;
        T w6 = w4 * (T(1.0) - t) + w5 * t;
        ctrlptsFront.emplace_back(p6 / w6);
        ctrlptsBack.emplace_back(p6 / w6);
        std::reverse(ctrlptsBack.begin(), ctrlptsBack.end());
        curvesOut
            = { BezierCurve<T>(ctrlptsFront,{c.m_weights[0],w1,w4,w6}, 0, 1, c.m_id), BezierCurve<T>(ctrlptsBack,{w6,w5,w3,c.m_weights[3]}, 0, 1, c.m_id) };
    }
    T tMid = t * (c.m_origTEnd - c.m_origTStart) + c.m_origTStart;
    curvesOut[0].m_origTStart = c.m_origTStart;
    curvesOut[0].m_origTEnd = tMid;
    curvesOut[0].m_origId = c.m_origId;
    curvesOut[1].m_origTStart = tMid;
    curvesOut[1].m_origTEnd = c.m_origTEnd;
    curvesOut[1].m_origId = c.m_origId;

    curvesOut[0].m_origLen = c.m_origLen;
    curvesOut[1].m_origLen = c.m_origLen;
    return true;
}

/**
 * @brief Splits a given curve at the specified t values while snapping the resulting
 *        curves endpoints to the precomputed curve(t) 2D point.
 *        Yields (number of splits + 1) curve pieces.
 *
 * @param curve curve to split
 * @param splitSet set of unique pairs <param t, precomputed curve(t)> as pair<T, Vec2<T>>
 * @param curvesOut resulting (number of splits + 1) curve pieces (output)
 * @return true if any splitting was performed
 * @return false else
 */
template <typename T>
bool split(const BezierCurve<T>& curve,
           const SetOfSplits& splitSet,
           vector<BezierCurve<T>>& curvesOut);

/**
 * @brief Split up all curves not well guardable
 *
 * @param curves curves to split (input and output)
 */
template <typename T>
void splitNonWellGuardables(vector<BezierCurve<T>>& curves);

/**
 * @brief Splits the given guarded bezier curves up if the current guarding quads are intersecting
 *        and recomputes the guards for the resulting curves. Curves occupying higher area will be
 *        split first.
 *
 * @param guardedCurves guarded curves to be split on intersections (input and output)
 * @return false if an unresolvable intersection was encountered and splitting was aborted
 * @return true else
 */
template <typename T>
bool splitIntersecting(vector<GuardedBezierCurve<T>>& guardedCurves);

} // namespace bzmsh

#endif
