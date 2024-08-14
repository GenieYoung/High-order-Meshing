#ifndef BZMSH_BEZIER_INJECTIVITY_CHECKER_HPP
#define BZMSH_BEZIER_INJECTIVITY_CHECKER_HPP

#include <queue>
#include <vector>

namespace bzmsh
{

/**
 * @brief Utility class to check validity of bezier triangles. 
 *
 */
template <typename T>
class BezierTriangleInjectivityChecker
{
  private:
    T min_area = 0.0;
    #ifdef ALL_EXACT
      int max_subdiv = 100; // number of subdivisions for conservative injectivity test
    #else
      int max_subdiv = 15;  // more does not make much sense in double precision
    #endif

    int current_level = 0;
    std::vector<T> v1, v2, w;
    int degreeB;
    int degreeC;
    int degreeRC;

    //(degree+1)*i - ((i*(i-1))/2) + j;
    std::vector<int> indB;
    std::vector<int> indC;
    std::vector<int> indRC;

    T b1(int i, int j)
    {
        return v1[indB[i] + j];
    }
    T b2(int i, int j)
    {
        return v2[indB[i] + j];
    }

    T xu(int i, int j)
    {
        return (b1(i + 1, j) - b1(i, j));
    } // constant factor does not matter: *3.0
    T xv(int i, int j)
    {
        return (b1(i, j + 1) - b1(i, j));
    } // constant factor does not matter: *3.0
    T yu(int i, int j)
    {
        return (b2(i + 1, j) - b2(i, j));
    } // constant factor does not matter: *3.0
    T yv(int i, int j)
    {
        return (b2(i, j + 1) - b2(i, j));
    } // constant factor does not matter: *3.0

    T P(const std::vector<T>& C, int i, int j, int n1, int n2, int n3);

    void bisect(const std::vector<T>& C,
                int degree,
                const std::vector<int>& ind,
                std::vector<T>& C0,
                std::vector<T>& C1);
    
    T facT(int n)
    {
        if (n == 0)
            return 1.0;
        if (n == 1)
            return 1.0;
        if (n == 2)
            return 2.0;
        if (n == 3)
            return 6.0;
        if (n == 4)
            return 24.0;
        if (n == 5)
            return 120.0;
        return fac(n - 1) * n;
    }

  public:
    BezierTriangleInjectivityChecker()
    {
    }

    // control points are expected in lexicographic index order, i.e b003, b012, b021, b030, b102,
    // b111, b120, b201, b210, b300.
    BezierTriangleInjectivityChecker(int _degree) : degreeB(_degree), degreeC(2 * (_degree - 1))
    {
        set_degree(_degree);
    }

    void set_degree(int degree);
    void set_max_subdiv(int subdiv)
    {
        max_subdiv = subdiv;
    }

    int is_definitely_injective_by_determinant(const std::vector<T>& _b1,
                                                const std::vector<T>& _b2,
                                                bool via_bisect = true);


    int is_definitely_injective_by_determinant(const std::vector<T>& _b1,
        const std::vector<T>& _b2,
        const std::vector<T>& _w,
        bool via_bisect = true);

    int level()
    {
        return current_level;
    }

    std::vector<T> derivate_u(const std::vector<T>& x)
    {
        std::vector<T> result;
        result.resize((degreeB) * (degreeB + 1) / 2.0);
        int k = 0;
        for (int i = 0; i <= degreeB - 1; ++i)
        {
            for (int j = 0; j <= degreeB - 1 - i; ++j)
            {
                result[k++] = x[indB[i + 1] + j] - x[indB[i] + j];
            }
        }
        return result;
    }
    std::vector<T> derivate_v(const std::vector<T>& x)
    {
        std::vector<T> result;
        result.resize((degreeB) * (degreeB + 1) / 2.0);
        int k = 0;
        for (int i = 0; i <= degreeB - 1; ++i)
        {
            for (int j = 0; j <= degreeB - 1 - i; ++j)
            {
                result[k++] = x[indB[i] + j + 1] - x[indB[i] + j];
            }
        }
        return result;
    }
    std::vector<T> multiply(const std::vector<T>& x1, const std::vector<T>& x2, int degree1_, int degree2_)
    {
        int degree1 = degree1_;
        int degree2 = degree2_;
        std::vector<T>result;
        int result_degree = degree1 + degree2;
        result.resize((result_degree + 1) * (result_degree + 2) / 2, T(0.0));

        std::vector<int>indC1(degree1 + 1), indC2(degree2 + 1), indR(result_degree + 1);
        int val = 0;
        int add = degree1 + 1;
        for (int i = 0; i <= degree1; ++i)
        {
            indC1[i] = val;
            val += add;
            add--;
        }
        val = 0;
        add = degree2 + 1;
        for (int i = 0; i <= degree2; ++i)
        {
            indC2[i] = val;
            val += add;
            add--;
        }
        val = 0;
        add = result_degree + 1;
        for (int i = 0; i <= result_degree; ++i)
        {
            indR[i] = val;
            val += add;
            add--;
        }

        for (int i = 0; i <= result_degree; i++)
            for (int j = 0; j <= result_degree - i; j++)
            {
                for (int ri = 0; ri <= degree2; ri++)
                    for (int rj = 0; rj <= degree2 - ri; rj++)
                    {
                        int si = i - ri;
                        int sj = j - rj;
                        int rk = degree2 - ri - rj;
                        int sk = degree1 - si - sj;
                        if (si < 0)
                            continue;
                        if (sj < 0)
                            continue;
                        if (si + sj > degree1)
                            continue;
                        if (result_degree - ri - rj - si - sj != result_degree - i - j)
                            continue;
                        T lijk = (facT(i) * facT(j) * facT(result_degree - i - j))
                            / (facT(ri) * facT(rj) * facT(rk) * facT(si) * facT(sj)
                                * facT(sk));
                        result[indR[i] + j] += lijk * (x2[indC2[ri] + rj] * x1[indC1[si] + sj]);
                    }
            }
        return result;
    }
    T facT(int degree, int i, int j)
    {
        if (i + j > degree || i > degree || j > degree)
        {
            return 0;
        }
        return facT(degree) / (facT(i) * facT(j) * facT(degree - i - j));
    }
};

} // namespace bzmsh

#include "util/bezier_injectivity_checker.cpp"

#endif
