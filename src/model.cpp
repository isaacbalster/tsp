/*
 *
 * tsp.cpp
 *
 * Traveling salesmen problem (TSP).
 *
 * Created by Isaac Balster on 26/01/2023.
 *
 */

#include "model.h"

std::vector<bool> getSubtour(const std::vector<std::vector<double>> & solution, const bool & silent,
                             const int & seed)
{
    int currentVertex = seed;
    bool closedLoop = false;

    auto nbVertices = static_cast<int>(solution.size());

    std::vector<bool> tour(nbVertices, false);

    /// We start the tour with currentVertex equal to the seed, which defaults to 0.
    tour[currentVertex] = true;

    if (!silent)
    {
        std::cout << "(integer) Subtour: " << seed << " ";
    }
    while (!closedLoop)
    {
        for (int j = 0; j < nbVertices; ++j)
        {
            if (solution[currentVertex][j] > TOLERANCE)
            {
                if (j == seed)
                {
                    closedLoop = true;
                }
                else
                {
                    currentVertex = j;
                    tour[currentVertex] = true;
                    if (!silent)
                    {
                        std::cout <<  currentVertex << " ";
                    }
                }
            }
        }
    }
    if (!silent)
    {
        std::cout << std::endl;
    }

    return tour;
}

void mtzModel(const ATSPDataC & data)
{
    int nbVertices = data._size;

    /// Passing pointers to gurobi environment and variables.
    GRBEnv *gurobiEnvironment = NULL;
    GRBVar **x = NULL;
    GRBVar *u = NULL;

    x = new GRBVar*[nbVertices];
    for (int i = 0; i < nbVertices; ++i)
    {
        x[i] = new GRBVar[nbVertices];
    }

    u = new GRBVar[nbVertices];

    try
    {
        gurobiEnvironment = new GRBEnv();
        GRBModel mtzModel = GRBModel(*gurobiEnvironment);

        /// Setting LazyConstraints parameter.
        mtzModel.set(GRB_IntParam_LazyConstraints, 1);

        /// Creating binary decision variables x_{i}_{j}.
        for (int i = 0; i < nbVertices; ++i)
        {
            for (int j = 0; j < nbVertices; ++j)
            {
                if (i != j)
                {
                    x[i][j] = mtzModel.addVar(0.0, 1.0, data._distances[i][j], 'B',
                                              "x_"+std::to_string(i)+"_"+std::to_string(j));
                }
                else
                    x[i][j] = mtzModel.addVar(0.0, 0.0, data._distances[i][j], 'I',
                                              "x_"+std::to_string(i)+"_"+std::to_string(j));
            }
        }

        /// Creating continuous decision variables u_{i}.
        for (int i = 0; i < nbVertices; ++i)
        {
            if (i == 0)
            {
                u[i] = mtzModel.addVar(0.0, 0.0, 0.0, 'C', "u_"+std::to_string(i));
            }
            else
            {
                u[i] = mtzModel.addVar(1.0, nbVertices-1, 0.0, 'C', "u_"+std::to_string(i));
            }
        }

        /// Degree constraints ($\sum_{j\in V:j\noteq i}x_{ij} = 1, \forall i\in V$ , i.e. one arc leaving)
        /// for each vertex $v\in V$.
        for (int i = 0; i < nbVertices; ++i)
        {
            GRBLinExpr expression = 0;
            for (int j = 0; j < nbVertices; ++j)
            {
                if (i != j)
                {
                    expression += x[i][j];
                }
            }
            mtzModel.addConstr(expression == 1, "degree_out_"+std::to_string(i));
        }

        /// Degree constraints ($\sum_{j\in V:j\noteq i}x_{ji} = 1, \forall i\in V$ , i.e. one arc entering)
        /// for each vertex $v\in V$.
        for (int i = 0; i < nbVertices; ++i)
        {
            GRBLinExpr expression = 0;
            for (int j = 0; j < nbVertices; ++j)
            {
                if (i != j)
                {
                    expression += x[j][i];
                }
            }
            mtzModel.addConstr(expression == 1, "degree_in_"+std::to_string(i));
        }

        for (int i = 0; i < nbVertices; ++i)
        {
            for (int j = 0; j < nbVertices; ++j)
            {
                if ((i != j) && (j != 0))
                {
                    GRBLinExpr expression = u[j] - u[i] - 1 + (nbVertices-1)*(1-x[i][j]);
                    mtzModel.addConstr(expression >= 0, "u_"+std::to_string(j)+"_"+std::to_string(i));
                }
            }
        }

        /// Writing the mathematical model in an *.lp file.
        mtzModel.write("mtz_tsp.lp");

        /// Optimize model
        mtzModel.optimize();

        // Extract solution
        if (mtzModel.get(GRB_IntAttr_SolCount) > 0)
        {
            /// Querying the optimal or best found solution.
            double **sol = new double*[nbVertices];
            for (int i = 0; i < nbVertices; ++i)
            {
                sol[i] = mtzModel.get(GRB_DoubleAttr_X, x[i], nbVertices);
            }

            int currentVertex = 0;
            int subtourLength = 1;
            bool closedLoop = false;

            std::cout << "Optimal tour: " << currentVertex << " ";
            while (!closedLoop)
            {
                for (int j = 0; j < nbVertices; ++j)
                {
                    if (sol[currentVertex][j] > TOLERANCE)
                    {
                        if (j == 0)
                        {
                            closedLoop = true;
                        }
                        else
                        {
                            currentVertex = j;
                            subtourLength += 1;
                            std::cout <<  currentVertex << " ";
                        }
                    }
                }
            }
            std::cout << std::endl;

            assert(subtourLength == nbVertices);

            for (int i = 0; i < nbVertices; i++)
                delete[] sol[i];
            delete[] sol;
        }

    } catch (GRBException e) {
        std::cout << "Error number: " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    } catch (...) {
        std::cout << "Error during optimization" << std::endl;
    }

    for (int i = 0; i < nbVertices; i++)
        delete[] x[i];
    delete[] x;
    delete gurobiEnvironment;
}

void lazyModel(const ATSPDataC & data, const std::string & lazyType, const bool & silent)
{
    int nbVertices = data._size;

    /// Passing pointers to gurobi environment and variables.
    GRBEnv *gurobiEnvironment = NULL;
    GRBVar **x = NULL;

    x = new GRBVar*[nbVertices];
    for (int i = 0; i < nbVertices; ++i)
    {
        x[i] = new GRBVar[nbVertices];
    }

    try
    {
        gurobiEnvironment = new GRBEnv();
        GRBModel lazyModel = GRBModel(*gurobiEnvironment);

        /// Setting LazyConstraints parameter.
        lazyModel.set(GRB_IntParam_LazyConstraints, 1);

        /// Creating binary decision variables x_{i}_{j}.
        for (int i = 0; i < nbVertices; ++i)
        {
            for (int j = 0; j < nbVertices; ++j)
            {
                if (i != j)
                {
                    x[i][j] = lazyModel.addVar(0.0, 1.0, data._distances[i][j], 'B',
                                              "x_"+std::to_string(i)+"_"+std::to_string(j));
                }
                else
                    x[i][j] = lazyModel.addVar(0.0, 0.0, data._distances[i][j], 'I',
                                              "x_"+std::to_string(i)+"_"+std::to_string(j));
            }
        }

        /// Degree constraints ($\sum_{j\in V:j\noteq i}x_{ij} = 1, \forall i\in V$ , i.e. one arc leaving)
        /// for each vertex $v\in V$.
        for (int i = 0; i < nbVertices; ++i)
        {
            GRBLinExpr expression = 0;
            for (int j = 0; j < nbVertices; ++j)
            {
                if (i != j)
                {
                    expression += x[i][j];
                }
            }
            lazyModel.addConstr(expression == 1, "degree_out_"+std::to_string(i));
        }

        /// Degree constraints ($\sum_{j\in V:j\noteq i}x_{ji} = 1, \forall i\in V$ , i.e. one arc entering)
        /// for each vertex $v\in V$.
        for (int i = 0; i < nbVertices; ++i)
        {
            GRBLinExpr expression = 0;
            for (int j = 0; j < nbVertices; ++j)
            {
                if (i != j)
                {
                    expression += x[j][i];
                }
            }
            lazyModel.addConstr(expression == 1, "degree_in_"+std::to_string(i));
        }

        /// Writing the mathematical model in an *.lp file.
        lazyModel.write("lazy_tsp.lp");

        /// Set callback function
        subtourElimination callback(x, nbVertices, lazyType, silent);
        lazyModel.setCallback(&callback);

        /// Optimize model
        lazyModel.optimize();

        // Extract solution
        if (lazyModel.get(GRB_IntAttr_SolCount) > 0)
        {
            /// Querying the optimal or best found solution.
            double **sol = new double*[nbVertices];
            for (int i = 0; i < nbVertices; ++i)
            {
                sol[i] = lazyModel.get(GRB_DoubleAttr_X, x[i], nbVertices);
            }

            int currentVertex = 0;
            int subtourLength = 1;
            bool closedLoop = false;

            std::cout << "Optimal tour: " << currentVertex << " ";
            while (!closedLoop)
            {
                for (int j = 0; j < nbVertices; ++j)
                {
                    if (sol[currentVertex][j] > TOLERANCE)
                    {
                        if (j == 0)
                        {
                            closedLoop = true;
                        }
                        else
                        {
                            currentVertex = j;
                            subtourLength += 1;
                            std::cout <<  currentVertex << " ";
                        }
                    }
                }
            }
            std::cout << std::endl;

            assert(subtourLength == nbVertices);

            for (int i = 0; i < nbVertices; i++)
                delete[] sol[i];
            delete[] sol;
        }

    } catch (GRBException e) {
        std::cout << "Error number: " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    } catch (...) {
        std::cout << "Error during optimization" << std::endl;
    }

    for (int i = 0; i < nbVertices; i++)
        delete[] x[i];
    delete[] x;
    delete gurobiEnvironment;
}
