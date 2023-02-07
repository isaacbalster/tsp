/*
 *
 * tsp.cpp
 *
 * Traveling salesmen problem (TSP).
 *
 * Created by Isaac Balster on 26/01/2023.
 *
 */

#ifndef _TSP_MODEL_H_
#define _TSP_MODEL_H_

#define TOLERANCE 1e-6

#include "data.h"

#include "hi_pr.hpp"
#include "gurobi_c++.h"

std::vector<bool> getSubtour(const std::vector<std::vector<double>> & solution, const bool & silent = true,
                             const int & seed = 0);

class subtourElimination: public GRBCallback
{
public:
    /// Members.
    int _n;
    std::string _lazyType;
    GRBVar** _vars;
    bool _silent;

    /// Constructor.
    subtourElimination(GRBVar **xVars, const int & nbVertices, std::string lazyType, const bool & silent = true) :
    _vars(xVars), _n(nbVertices), _lazyType(std::move(lazyType)), _silent(silent)
    {
        if ((_lazyType != "SEC") && (_lazyType != "CUT"))
        {
            std::cout << "Lazy constraints type not known. Constraints of type "
                      << "$\\sum_{i\\in N\\setminus S}\\sum_{j\\in S}x_{ij}\\geq 1$ will be added "
                      << "to ensure feasibility." << std::endl;
        }
        if ((_lazyType == "SEC"))
        {
            std::cout << "SEC lazy constraints. Constraints of type "
                      << "$\\sum_{i\\in S}\\sum_{j\\in S:i\\neq j}x_{ij}\\leq |S|-1$ are added "
                      << "to ensure feasibility." << std::endl;
        }
        if ((_lazyType == "CUT"))
        {
            std::cout << "CUT lazy constraints. Constraints of type "
                      << "$\\sum_{i\\in N\\setminus S}\\sum_{j\\in S}x_{ij}\\geq 1$ are added "
                      << "to ensure feasibility." << std::endl;
        }
    }

    /// Destructor.
    ~subtourElimination() override = default;

protected:
    void callback() override
    {
        try
        {
            if (where == GRB_CB_MIPNODE)
            {
                if (!_silent)
                {
                    std::cout << "*****  NEW ROUND OF FRACTIONAL CUT SEPARATION  *****" << std::endl;
                }

                double maxFlow;
                long *distance = new long[_n];

                for (int i = 0; i < _n; ++i)
                {
                    distance[i] = _n;
                }

                double **x = new double*[_n];
                for (int i = 0; i < _n; ++i)
                {
                    x[i] = getSolution(_vars[i], _n);
                }

                int nbNodesSwept = 0;
                int subtourLength = 0;
                int complementSubtourLength = _n;
                std::vector<bool> sweptNodes(_n, false);
                std::vector<bool> subtour(_n, false);
                std::vector<bool> complementSubtour(_n, true);
                sweptNodes[0] = true;

                int sinkNode = 1;

                /// Computation of the minimum cut separating 0 (source) of sinkNode (sink).
                while (nbNodesSwept < _n - 1)
                {
                    directed_min_cut(x,_n,0,sinkNode,maxFlow,distance);
                    sweptNodes[sinkNode] = true;
                    nbNodesSwept += 1;

                    if (maxFlow < 1 - TOLERANCE)
                    {
                        subtour[sinkNode] = true;
                        subtourLength = 0;

                        for (int i = 1; i < _n; ++i)
                        {
                            if (distance[i] < _n)
                            {
                                subtour[i] = true;
                                complementSubtour[i] = false;
                                sweptNodes[i] = true;
                                subtourLength += 1;
                                complementSubtourLength -= 1;
                                if (i != sinkNode)
                                {
                                    nbNodesSwept += 1;
                                }
                            }
                        }

                        if (subtourLength > 1)
                        {
                            if (!_silent)
                            {
                                //std::cout << "sinkNode is "<< sinkNode << " with maxFlow equal to " << maxFlow
                                //          << " and distances[" << _n << "] is equal " << "to [ ";
                                //for (int j = 0; j < _n; j++)
                                //{
                                //    std::cout << distance[j] << " ";
                                //}
                                //std::cout << "]. ";
                                std::cout << "(fractional) Subtour: ";
                                for (int j = 1; j < _n; j++)
                                {
                                    if (subtour[j])
                                        std::cout << j << " ";
                                }
                                std::cout << std::endl;
                            }

                            /// Adding subtour elimination constraint.
                            GRBLinExpr expression = 0;
                            double leftHandSide = 0.0f;

                            for (int i = 0; i < _n; i++)
                            {
                                for (int j = 0; j < _n; j++)
                                {
                                    if ((!subtour[i]) && (subtour[j]))
                                    {
                                        expression += _vars[i][j];
                                        if (x[i][j] > TOLERANCE)
                                        {
                                            leftHandSide += x[i][j];
                                        }
                                    }
                                }
                            }
                            if (1 - leftHandSide > TOLERANCE)
                            {
                                if (!_silent)
                                {
                                    std::cout << "User cut : ";
                                    for (int i = 0; i < _n; i++)
                                    {
                                        for (int j = 0; j < _n; j++)
                                        {
                                            if ((!subtour[i]) && (subtour[j]))
                                            {
                                                std::cout << "x[" << i << "][" << j << "]";
                                                if (x[i][j] > TOLERANCE)
                                                {
                                                    std::cout << "(" << x[i][j] << ") + ";
                                                }
                                                else
                                                {
                                                    std::cout << " + ";
                                                }
                                            }
                                        }
                                    }
                                    std::cout << " >= 1. lhs = " << leftHandSide << ", rhs = " << 1
                                              << ", violation = " <<  1 - leftHandSide << "." << std::endl;
                                }

                                addCut(expression >= 1);
                            }
                        }

                        std::fill(subtour.begin(), subtour.end(), false);
                        subtourLength = 0;
                    }
                    for (int i = 1; i < _n; ++i)
                    {
                        if (!sweptNodes[i])
                        {
                            sinkNode = i;
                            for (int j = 0; j < _n; ++j)
                            {
                                distance[j] = _n;
                            }
                            break;
                        }
                    }
                }

                if ((complementSubtourLength > 1) && (complementSubtourLength < _n))
                {
                    if (!_silent)
                    {
                        std::cout << "(fractional) Subtour: ";
                        for (int j = 0; j < _n; j++)
                        {
                            if (complementSubtour[j])
                                std::cout << j << " ";
                        }
                        std::cout << std::endl;
                    }

                    /// Adding subtour elimination constraint.
                    GRBLinExpr expression = 0;
                    double leftHandSide = 0.0f;

                    for (int i = 0; i < _n; i++)
                    {
                        for (int j = 0; j < _n; j++)
                        {
                            if ((!complementSubtour[i]) && (complementSubtour[j]))
                            {
                                expression += _vars[i][j];
                                if (x[i][j] > TOLERANCE)
                                {
                                    leftHandSide += x[i][j];
                                }
                            }
                        }
                    }
                    if (1 - leftHandSide > TOLERANCE)
                    {
                        if (!_silent)
                        {
                            std::cout << "User cut : ";
                            for (int i = 0; i < _n; i++)
                            {
                                for (int j = 0; j < _n; j++)
                                {
                                    if ((!complementSubtour[i]) && (complementSubtour[j]))
                                    {
                                        std::cout << "x[" << i << "][" << j << "]";
                                        if (x[i][j] > TOLERANCE)
                                        {
                                            std::cout << "(" << x[i][j] << ") + ";
                                        }
                                        else
                                        {
                                            std::cout << " + ";
                                        }
                                    }
                                }
                            }
                            std::cout << " >= 1. lhs = " << leftHandSide << ", rhs = " << 1
                                      << ", violation = " <<  1 - leftHandSide << "." << std::endl;
                        }

                        addCut(expression >= 1);
                    }
                }
            }

            if (where == GRB_CB_MIPSOL)
            {
                if (!_silent)
                {
                    std::cout << "*****  NEW ROUND OF INTEGER CUT SEPARATION  *****" << std::endl;
                }

                /// Found an integer feasible solution, query it.
                std::vector<std::vector<double>> solution(_n, std::vector<double>(_n, 0.0f));
                for (int i = 0; i < _n; ++i)
                {
                    for (int j = 0; j < _n; ++j)
                    {
                        solution[i][j] = getSolution(_vars[i][j]);
                    }
                }

                /// Then check if it visits every node.
                std::vector<bool> sweptNodes(_n, false);
                int nbSweptNodes = 0;
                int seed = 0;

                while (nbSweptNodes < _n)
                {
                    int subtourLength = 0;
                    std::vector<bool> subtour(getSubtour(solution, _silent, seed));

                    for (int i = 0; i < _n; i++)
                    {
                        if (subtour[i])
                        {
                            subtourLength += 1;
                            sweptNodes[i] = true;
                            nbSweptNodes += 1;
                        }
                    }

                    /// If not, we add a lazy constraint.
                    if (subtourLength < _n)
                    {
                        /// Adding subtour elimination constraints.
                        GRBLinExpr expression = 0;
                        double leftHandSide = 0.0f;

                        if (!_silent)
                        {
                            std::cout << "Lazy constraint : ";
                        }
                        if (_lazyType == "SEC")
                        {
                            for (int i = 0; i < _n; i++)
                            {
                                for (int j = 0; j < _n; j++)
                                {
                                    if ((i != j) && (subtour[i]) && (subtour[j]))
                                    {
                                        expression += _vars[i][j];
                                        if (!_silent)
                                        {
                                            std::cout << "x[" << i << "][" << j << "]";
                                            if (solution[i][j] > TOLERANCE)
                                            {
                                                std::cout << "(" << solution[i][j] << ") + ";
                                                leftHandSide += solution[i][j];
                                            }
                                            else
                                            {
                                                std::cout << " + ";
                                            }
                                        }
                                    }
                                }
                            }
                            addLazy(expression <= subtourLength - 1);
                            if (!_silent)
                            {
                                std::cout << " <= " << subtourLength - 1 << ". lhs = " << leftHandSide << ", rhs = "
                                          << subtourLength - 1 << ", violation = " << leftHandSide - (subtourLength - 1)
                                          << "." << std::endl;
                            }
                        }
                        else /// if _lazyType == "CUT" or user inputs a typo in lazyType parameter.
                        {
                            for (int i = 0; i < _n; i++)
                            {
                                for (int j = 0; j < _n; j++)
                                {
                                    if ((!subtour[i]) && (subtour[j]))
                                    {
                                        expression += _vars[i][j];
                                        if (!_silent)
                                        {
                                            std::cout << "x[" << i << "][" << j << "]";
                                            if (solution[i][j] > TOLERANCE)
                                            {
                                                std::cout << "(" << solution[i][j] << ") + ";
                                                leftHandSide += solution[i][j];
                                            }
                                            else
                                            {
                                                std::cout << " + ";
                                            }
                                        }
                                    }
                                }
                            }
                            addLazy(expression >= 1);
                            if (!_silent)
                            {
                                std::cout << " >= 1. lhs = " << leftHandSide << ", rhs = " << 1
                                          << ", violation = " <<  1 - leftHandSide << "." << std::endl;
                            }
                        }
                    }
                    while (sweptNodes[seed] && seed < _n)
                    {
                        ++seed;
                    }
                }
            }
        } catch (GRBException e)
        {
            std::cout << "Error number: " << e.getErrorCode() << std::endl;
            std::cout << e.getMessage() << std::endl;
        } catch (...)
        {
            std::cout << "Error during callback" << std::endl;
        }
    }
};

void mtzModel(const ATSPDataC & data);

void lazyModel(const ATSPDataC & data, const std::string & lazyType, const bool & silent = true);

#endif /// _TSP_MODEL_H_