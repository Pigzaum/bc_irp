////////////////////////////////////////////////////////////////////////////////
/*
 * File: init_grb_model.cpp
 * Author: Guilherme O. Chagas
 *
 * @brief Gurobi's model initialization helper functions definitions.
 *
 * @acknowledgment Special thanks to Ph.D. Cleder Marcos Schenekemberg.
 * 
 * (I'm sorry for my bad english xD)
 *
 * Created on October 22, 2020, 10:44 PM.
 * 
 * References:
 */
////////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <sstream>

#include "../include/ext/loguru/loguru.hpp"

#include "../include/init_grb_model.hpp"


namespace
{

void GetSubsetsSubtour(int i,
                       int num_vertex,
                       int &count_subset,
                       std::vector<int>& Set,
                       std::vector<std::vector<int>>& S)
{
    /* subsets elimination subtour S \in V (S != {}) */

    if (i > num_vertex)
    {
        for (int j = 1; j <= num_vertex; j++)
        {
            if (Set[j] == 1)
            {
                S[count_subset].push_back(j);
            }
        }
        count_subset++;
    }else
    {
        Set[i] = 1;
        GetSubsetsSubtour(i + 1, num_vertex, count_subset, Set, S);
        Set[i] = 0;
        GetSubsetsSubtour(i + 1, num_vertex, count_subset, Set, S);
    }
}

} // anonymous namespace


void init::inventoryLevelVariables(
    GRBModel& model,
    std::vector<std::vector<GRBVar>>& I,
    const std::shared_ptr<const Instance>& pInst)
{
    DRAW_LOG_F(INFO, "\tinitializing inventory level variables (I_it)");

    I.reserve(pInst->getNbVertices());
    for (auto i = 0; i < pInst->getNbVertices(); ++i)
    {
        I.push_back(std::vector<GRBVar>());
        I[i].reserve(pInst->getT() + 1);
        for (auto t = 0; t <= pInst->getT(); ++t)
        {
            std::ostringstream oss;
            oss << "I_" << i << "_" << t;
            I[i].push_back(model.addVar(0,
                                        i == 0 ? GRB_INFINITY : pInst->getUi(i),
                                        pInst->get_hi(i),
                                        GRB_CONTINUOUS,
                                        oss.str()));
        }
    }
}

void init::quantityVariables(GRBModel& model,
                             std::vector<std::vector<GRBVar>>& q,
                             const std::shared_ptr<const Instance>& pInst)
{
    DRAW_LOG_F(INFO, "\tinitializing quantity delivered variables (q_it)");

    q.reserve(pInst->getNbVertices());
    for (auto i = 0; i < pInst->getNbVertices(); ++i)
    {
        q.push_back(std::vector<GRBVar>(pInst->getT() + 1));

        if (i == 0) continue;

        for (auto t = 0; t < pInst->getT(); ++t)
        {
            std::ostringstream oss;
            oss << "q_" << i << "_" << t;
            q[i][t] = model.addVar(0,
                                   GRB_INFINITY,
                                   0,
                                   GRB_CONTINUOUS,
                                   oss.str());
        }
    }
}


void init::visitationVariables(GRBModel& model,
                               std::vector<std::vector<GRBVar>>& y,
                               const std::shared_ptr<const Instance>& pInst)
{
    DRAW_LOG_F(INFO, "\tinitializing visitation variables (y_it)");

    y.reserve(pInst->getNbVertices());
    for (auto i = 0; i < pInst->getNbVertices(); ++i)
    {
        y.push_back(std::vector<GRBVar>(pInst->getT() + 1));
        for (auto t = 0; t < pInst->getT(); ++t)
        {
            std::ostringstream oss;
            oss << "y_" << i << "_" << t;
            y[i][t] = model.addVar(0, 1, 0, GRB_BINARY, oss.str());
        }
    }
}


void init::routingVariables(GRBModel& model,
                            std::vector<std::vector<std::vector<GRBVar>>>& x,
                            const std::shared_ptr<const Instance>& pInst)
{
    DRAW_LOG_F(INFO, "\tinitializing routing variables (x_ijt)");

    x.reserve(pInst->getNbVertices());
    for (auto i = 0; i < pInst->getNbVertices(); ++i)
    {
        x.push_back(std::vector<std::vector<GRBVar>>());
        x[i].reserve(pInst->getNbVertices());
        for (auto j = 0; j < pInst->getNbVertices(); ++j)
        {
            x[i].push_back(std::vector<GRBVar>(pInst->getT() + 1));

            if (j <= i) continue;

            for (auto t = 0; t < pInst->getT(); ++t)
            {
                std::ostringstream oss;
                oss << "x_" << i << "_" << j << "_" << t;
                x[i][j][t] = model.addVar(0,
                                          (i == 0 ? 2 : 1),
                                          pInst->get_cij(i, j),
                                          GRB_INTEGER,
                                          oss.str());
            }
        }
    }
}


void init::inventoryDefDepotConstrs(
    GRBModel& model,
    std::vector<GRBConstr>& constrs,
    const std::vector<std::vector<GRBVar>>& I,
    const std::vector<std::vector<GRBVar>>& q,
    const std::shared_ptr<const Instance>& pInst)
{
    DRAW_LOG_F(INFO, "\tinitializing inventory definition depot constraints");

    constrs.push_back(model.addConstr(I[0][0] == pInst->getIi0(0), "1C_0"));
    for (auto t = 1; t <= pInst->getT(); ++t)
    {
        GRBLinExpr e = I[0][t - 1] + pInst->get_rit(0, t - 1);
        for (auto i = 1; i < pInst->getNbVertices(); ++i) // skip depot
        {
            e -= q[i][t - 1];
        }

        std::ostringstream oss;
        oss << "1C_" << t;
        constrs.push_back(model.addConstr(I[0][t] == e, oss.str()));
    }
}


void init::stockOutDepotConstrs(GRBModel& model,
                                std::vector<GRBConstr>& constrs,
                                const std::vector<std::vector<GRBVar>>& I,
                                const std::vector<std::vector<GRBVar>>& q,
                                const std::shared_ptr<const Instance>& pInst)
{
    DRAW_LOG_F(INFO, "\tinitializing stockout depot contraints");

    for (auto t = 0; t < pInst->getT(); ++t)
    {
        GRBLinExpr e = 0;
        for (int i = 1; i < pInst->getNbVertices(); i++)
        {
            e += q[i][t];
        }

        std::ostringstream oss;
        oss << "2C" << t;
        constrs.push_back(model.addConstr(I[0][t] >= e, oss.str()));
    }
}


void init::inventoryDefCustomersConstrs(
    GRBModel& model,
    std::vector<GRBConstr>& constrs,
    const std::vector<std::vector<GRBVar>>& I,
    const std::vector<std::vector<GRBVar>>& q,
    const std::shared_ptr<const Instance>& pInst)
{
    DRAW_LOG_F(INFO, "\tinitializing inventory definition customers contraints");

    for (auto i = 1; i < pInst->getNbVertices(); ++i)
    {
        std::ostringstream oss;
        oss << "3C_" << i;
        constrs.push_back(
            model.addConstr(I[i][0] == pInst->getIi0(i), oss.str()));

        for (auto t = 1; t <= pInst->getT(); ++t)
        {
            GRBLinExpr e = I[i][t - 1] + q[i][t - 1] - pInst->get_rit(i, t - 1);
            oss.clear();
            oss.str("");
            oss << "3C_" << i << "_" << t;
            constrs.push_back(model.addConstr(I[i][t] == e, oss.str()));
        }
    }
}


void init::capacityConstrs(GRBModel& model,
                           std::vector<GRBConstr>& constrs,
                           const std::vector<std::vector<GRBVar>>& q,
                           const std::shared_ptr<const Instance>& pInst)
{
    DRAW_LOG_F(INFO, "\tinitializing capacity contraints");

    for (auto t = 0; t < pInst->getT(); ++t)
    {
        GRBLinExpr e = 0;
        for (auto i = 1; i < pInst->getNbVertices(); ++i)
        {
            e += q[i][t];
        }

        std::ostringstream oss;
        oss << "4C_" << t;
        constrs.push_back(model.addConstr(e <= pInst->getC(), oss.str()));
    }
}


void init::inventoryLevelConstrs(
    GRBModel& model,
    std::vector<GRBConstr>& constrs,
    const std::vector<std::vector<GRBVar>>& I,
    const std::vector<std::vector<GRBVar>>& q,
    const std::shared_ptr<const Instance>& pInst)
{
    DRAW_LOG_F(INFO, "\tinitializing inventory level constraints");

    for (int t = 0; t < pInst->getT(); ++t)
    {
        GRBLinExpr e = 0;
        for (int i = 1; i < pInst->getNbVertices(); ++i)
        {
            e += q[i][t];
        }

        std::ostringstream oss;
        oss << "5C_" << t;
        constrs.push_back(model.addConstr(I[0][t] >= e, oss.str()));
    }
}


void init::mlQuantityCapacityConstrs(
    GRBModel& model,
    std::vector<GRBConstr>& constrs,
    const std::vector<std::vector<GRBVar>>& I,
    const std::vector<std::vector<GRBVar>>& q,
    const std::shared_ptr<const Instance>& pInst)
{
    DRAW_LOG_F(INFO, "\tML quantity capacity constraints");

    for (auto i = 1; i < pInst->getNbVertices(); ++i)
    {
        for (auto t = 0; t < pInst->getT(); ++t)
        {
            std::ostringstream oss;
            oss << "6C_ML" << i << "_" << t;
            constrs.push_back(
                model.addConstr(q[i][t] <= pInst->getUi(i) - I[i][t],
                oss.str()));
        }
    }
}


void init::ouQuantityCapacityConstrs(
    GRBModel& model,
    std::vector<GRBConstr>& constrs,
    const std::vector<std::vector<GRBVar>>& I,
    const std::vector<std::vector<GRBVar>>& y,
    const std::vector<std::vector<GRBVar>>& q,
    const std::shared_ptr<const Instance>& pInst)
{
    DRAW_LOG_F(INFO, "\tOU quantity capacity constraints");

    for (auto i = 1; i < pInst->getNbVertices(); ++i)
    {
        for (auto t = 0; t < pInst->getT(); ++t)
        {
            std::ostringstream oss;
            oss << "6C_OU" << i << "_" << t;
            constrs.push_back(model.addConstr(q[i][t] >=
                              pInst->getUi(i) * y[i][t] - I[i][t], oss.str()));
        }
    }
}


void init::quantitiesRoutingConstraint(
    GRBModel& model,
    std::vector<GRBConstr>& constrs,
    const std::vector<std::vector<GRBVar>>& y,
    const std::vector<std::vector<GRBVar>>& q,
    const std::shared_ptr<const Instance>& pInst)
{
    DRAW_LOG_F(INFO, "\tinitializing quantities routing constraints");

    for (auto i = 1; i < pInst->getNbVertices(); ++i)
    {
        for (auto t = 0; t < pInst->getT(); ++t)
        {
            std::ostringstream oss;
            oss << "7C" << i << "_" << t;
            constrs.push_back(
                model.addConstr(q[i][t] <= pInst->getUi(i) * y[i][t],
                oss.str()));
        }
    }
}


void init::capacityVehicleConstraint(
    GRBModel& model,
    std::vector<GRBConstr>& constrs,
    const std::vector<std::vector<GRBVar>>& y,
    const std::vector<std::vector<GRBVar>>& q,
    const std::shared_ptr<const Instance>& pInst)
{
    DRAW_LOG_F(INFO, "\tinitializing capacity vehicles constraints");

    for (auto t = 0; t < pInst->getT(); ++t)
    {
        GRBLinExpr e = 0;
        for (auto i = 1; i < pInst->getNbVertices(); ++i)
        {
            e += q[i][t];
        }

        std::ostringstream oss;
        oss << "8C_" << t;
        constrs.push_back(
            model.addConstr(e <= pInst->getC() * y[0][t], oss.str()));
    }
}


void init::degreeConstrs(GRBModel& model,
                         std::vector<GRBConstr>& constrs,
                         const std::vector<std::vector<GRBVar>>& y,
                         const std::vector<std::vector<std::vector<GRBVar>>>& x,
                         const std::shared_ptr<const Instance>& pInst)
{
    DRAW_LOG_F(INFO, "\tinitializing degree constraints");

    for (auto t = 0; t < pInst->getT(); ++t)
    {
        for (auto i = 0; i < pInst->getNbVertices(); ++i)
        {
            GRBLinExpr e1 = 0;
            for (auto j = i + 1; j < pInst->getNbVertices(); ++j)
            {
                e1 += x[i][j][t];
            }

            GRBLinExpr e2 = 0;
            for (auto j = 0; j < i; ++j)
            {
                e2 += x[j][i][t];
            }

            std::ostringstream oss;
            oss << "9C_" << i << "_" << t;
            constrs.push_back(
                model.addConstr(e1 + e2 == 2 * y[i][t], oss.str()));
        }
    }
}


void init::subtourEliminationConstrs(
    GRBModel& model,
    std::vector<GRBConstr>& constrs,
    const std::vector<std::vector<GRBVar>>& y,
    const std::vector<std::vector<std::vector<GRBVar>>>& x,
    const std::shared_ptr<const Instance>& pInst)
{
    DRAW_LOG_F(INFO, "\tinitializing subtour elimination constraints");

    /* subsets */
    int count_subset = 0;

    std::vector<int> Set(pInst->getNbVertices());
    int number_of_subsets = std::pow(2, pInst->getNbVertices() - 1) - 1;

    /* V' = {1, ..., num_vertex}. S contains the subsets of V (S != {}) */
    std::vector<std::vector<int>> S(number_of_subsets);
    GetSubsetsSubtour(1, pInst->getNbVertices() - 1, count_subset, Set, S);

    /* for each subset */
    for (int subset = 0; subset < number_of_subsets; subset++)
    {
        if (static_cast<int>(S[subset].size()) > 1)
        {
            for (auto t = 0; t < pInst->getT(); ++t)
            {
                for (size_t m = 0; m < S[subset].size(); m++)
                {
                    GRBLinExpr lhs = 0;
                    GRBLinExpr rhs = 0;

                    /* left expression */
                    for (size_t i = 0; i < S[subset].size(); i++)
                    {
                        for (size_t j = 0; j < S[subset].size(); j++)
                        {
                            if (S[subset][i] < S[subset][j])
                            {
                                lhs += x[S[subset][i]][S[subset][j]][t];
                            }
                        }
                    }

                    /* right expression */
                    for (size_t i = 0; i < S[subset].size(); i++)
                    {
                        rhs += y[S[subset][i]][t];
                    }

                    std::ostringstream oss;
                    oss << "10C" << subset+1 << "_" << t << "_" << S[subset][m];
                    rhs -= y[S[subset][m]][t];
                    constrs.push_back(model.addConstr(lhs <= rhs, oss.str()));
                }
            }
        }
    }
}