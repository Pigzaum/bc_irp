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

#include <bitset>
#include <cmath>
#include <sstream>

#include "../include/ext/loguru/loguru.hpp"

#include "../include/init_grb_model.hpp"

/////////////////////////////// Helper functions ///////////////////////////////

namespace
{

std::vector<int> findSetS(const int idx, const int nbVertices)
{
    const int MAX_N = 30;
    CHECK_F(nbVertices <= MAX_N);

    std::vector<int> inS;
    inS.reserve(nbVertices);

    std::bitset<MAX_N> myBitset(idx);
    for (int i = 1; i < nbVertices; ++i)
    {
        if (myBitset[i - 1])
        {
            inS.push_back(i);
        }
    }

    return inS;
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
                             std::vector<std::vector<std::vector<GRBVar>>>& q,
                             const std::shared_ptr<const Instance>& pInst)
{
    DRAW_LOG_F(INFO, "\tinitializing quantity delivered variables (q_it)");

    q.reserve(pInst->getNbVertices());
    for (auto i = 0; i < pInst->getNbVertices(); ++i)
    {
        q.push_back(std::vector<std::vector<GRBVar>>());
        q[i].reserve(pInst->getK());

        for (auto k = 0; k < pInst->getK(); ++k)
        {
            q[i].push_back(std::vector<GRBVar>(pInst->getT() + 1));

            if (i == 0) continue;

            for (auto t = 0; t < pInst->getT(); ++t)
            {
                std::ostringstream oss;
                oss << "q_" << i << "_" << k << "_" << t;
                q[i][k][t] = model.addVar(0,
                                    GRB_INFINITY,
                                    0,
                                    GRB_CONTINUOUS,
                                    oss.str());
            }
        }
    }
}


void init::visitationVariables(GRBModel& model,
                               std::vector<std::vector<std::vector<GRBVar>>>& y,
                               const std::shared_ptr<const Instance>& pInst)
{
    DRAW_LOG_F(INFO, "\tinitializing visitation variables (y_it)");

    y.reserve(pInst->getNbVertices());
    for (auto i = 0; i < pInst->getNbVertices(); ++i)
    {
        y.push_back(std::vector<std::vector<GRBVar>>());
        y[i].reserve(pInst->getK());
        for (auto k = 0; k < pInst->getK(); ++k)
        {
            y[i].push_back(std::vector<GRBVar>(pInst->getT() + 1));
            for (auto t = 0; t < pInst->getT(); ++t)
            {
                std::ostringstream oss;
                oss << "y_" << i << "_" << k << "_" << t;
                y[i][k][t] = model.addVar(0, 1, 0, GRB_BINARY, oss.str());
            }
        }
    }
}


void init::routingVariables(
    GRBModel& model,
    std::vector<std::vector<std::vector<std::vector<GRBVar>>>>& x,
    const std::shared_ptr<const Instance>& pInst)
{
    DRAW_LOG_F(INFO, "\tinitializing routing variables (x_ijt)");

    x.reserve(pInst->getNbVertices());
    for (auto i = 0; i < pInst->getNbVertices(); ++i)
    {
        x.push_back(std::vector<std::vector<std::vector<GRBVar>>>());
        x[i].reserve(pInst->getNbVertices());
        for (auto j = 0; j < pInst->getNbVertices(); ++j)
        {
            x[i].push_back(std::vector<std::vector<GRBVar>>());
            x[i][j].reserve(pInst->getK());

            for (auto k = 0; k < pInst->getK(); ++k)
            {
                x[i][j].push_back(std::vector<GRBVar>(pInst->getT() + 1));
                if (j <= i) continue;

                for (auto t = 0; t < pInst->getT(); ++t)
                {
                    std::ostringstream oss;
                    oss << "x_" << i << "_" << j << "_" << k << "_" << t;
                    x[i][j][k][t] = model.addVar(0,
                                                 (i == 0 ? 2 : 1),
                                                 pInst->get_cij(i, j),
                                                 GRB_INTEGER,
                                                 oss.str());
                }
            }
        }
    }
}


// ok
void init::inventoryDefDepotConstrs(
    GRBModel& model,
    std::vector<GRBConstr>& constrs,
    const std::vector<std::vector<GRBVar>>& I,
    const std::vector<std::vector<std::vector<GRBVar>>>& q,
    const std::shared_ptr<const Instance>& pInst)
{
    DRAW_LOG_F(INFO, "\tinitializing inventory definition depot constraints");

    constrs.push_back(model.addConstr(I[0][0] == pInst->getIi0(0), "1C_0"));
    for (auto t = 1; t <= pInst->getT(); ++t)
    {
        GRBLinExpr e = I[0][t - 1] + pInst->get_rit(0, t - 1);
        for (auto i = 1; i < pInst->getNbVertices(); ++i) // skip depot
        {
            for (auto k = 0; k < pInst->getK(); ++k)
            {
                e -= q[i][k][t - 1];
            }
        }

        std::ostringstream oss;
        oss << "1C_" << t;
        constrs.push_back(model.addConstr(I[0][t] == e, oss.str()));
    }
}


// ok
void init::stockOutDepotConstrs(
    GRBModel& model,
    std::vector<GRBConstr>& constrs,
    const std::vector<std::vector<GRBVar>>& I,
    const std::vector<std::vector<std::vector<GRBVar>>>& q,
    const std::shared_ptr<const Instance>& pInst)
{
    DRAW_LOG_F(INFO, "\tinitializing stockout depot contraints");

    for (auto t = 0; t < pInst->getT(); ++t)
    {
        GRBLinExpr e = 0;
        for (int i = 1; i < pInst->getNbVertices(); i++)
        {
            for (auto k = 0; k < pInst->getK(); ++k)
            {
                e += q[i][k][t];
            }
        }

        std::ostringstream oss;
        oss << "2C_" << t;
        constrs.push_back(model.addConstr(I[0][t] >= e, oss.str()));
    }
}


// ok
void init::inventoryDefCustomersConstrs(
    GRBModel& model,
    std::vector<GRBConstr>& constrs,
    const std::vector<std::vector<GRBVar>>& I,
    const std::vector<std::vector<std::vector<GRBVar>>>& q,
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
            GRBLinExpr e = 0;
            for (auto k = 0; k < pInst->getK(); ++k)
            {
                e += q[i][k][t - 1];
            }

            e += I[i][t - 1] - pInst->get_rit(i, t - 1);

            oss.clear();
            oss.str("");
            oss << "3C_" << i << "_" << t;
            constrs.push_back(model.addConstr(I[i][t] == e, oss.str()));
        }
    }
}


// ok
void init::capacityConstrs(
    GRBModel& model,
    std::vector<GRBConstr>& constrs,
    const std::vector<std::vector<std::vector<GRBVar>>>& q,
    const std::shared_ptr<const Instance>& pInst)
{
    DRAW_LOG_F(INFO, "\tinitializing capacity contraints");

    for (auto t = 0; t < pInst->getT(); ++t)
    {
        GRBLinExpr e = 0;
        for (auto i = 1; i < pInst->getNbVertices(); ++i)
        {
            for (auto k = 0; k < pInst->getK(); ++k)
            {
                e += q[i][k][t];
            }
        }

        std::ostringstream oss;
        oss << "4C_" << t;
        constrs.push_back(model.addConstr(e <= pInst->getC(), oss.str()));
    }
}


// ok
void init::inventoryLevelConstrs(
    GRBModel& model,
    std::vector<GRBConstr>& constrs,
    const std::vector<std::vector<GRBVar>>& I,
    const std::vector<std::vector<std::vector<GRBVar>>>& q,
    const std::shared_ptr<const Instance>& pInst)
{
    DRAW_LOG_F(INFO, "\tinitializing inventory level constraints");

    for (int t = 0; t < pInst->getT(); ++t)
    {
        GRBLinExpr e = 0;
        for (int i = 1; i < pInst->getNbVertices(); ++i)
        {
            for (auto k = 0; k < pInst->getK(); ++k)
            {
                e += q[i][k][t];
            }
        }

        std::ostringstream oss;
        oss << "5C_" << t;
        constrs.push_back(model.addConstr(I[0][t] >= e, oss.str()));
    }
}


// ok
void init::mlQuantityCapacityConstrs(
    GRBModel& model,
    std::vector<GRBConstr>& constrs,
    const std::vector<std::vector<GRBVar>>& I,
    const std::vector<std::vector<std::vector<GRBVar>>>& q,
    const std::shared_ptr<const Instance>& pInst)
{
    DRAW_LOG_F(INFO, "\tML quantity capacity constraints");

    for (auto i = 1; i < pInst->getNbVertices(); ++i)
    {
        for (auto t = 0; t < pInst->getT(); ++t)
        {
            GRBLinExpr e = 0;
            for (auto k = 0; k < pInst->getK(); ++k)
            {
                e += q[i][k][t];
            }

            std::ostringstream oss;
            oss << "6C_ML" << i << "_" << t;
            constrs.push_back(
                model.addConstr(e <= pInst->getUi(i) - I[i][t], oss.str()));
        }
    }
}


// ok
void init::ouQuantityCapacityConstrs(
    GRBModel& model,
    std::vector<GRBConstr>& constrs,
    const std::vector<std::vector<GRBVar>>& I,
    const std::vector<std::vector<std::vector<GRBVar>>>& y,
    const std::vector<std::vector<std::vector<GRBVar>>>& q,
    const std::shared_ptr<const Instance>& pInst)
{
    DRAW_LOG_F(INFO, "\tOU quantity capacity constraints");

    for (auto i = 1; i < pInst->getNbVertices(); ++i)
    {
        for (auto t = 0; t < pInst->getT(); ++t)
        {
            GRBLinExpr e1 = 0;
            GRBLinExpr e2 = 0;
            for (int k = 0; k < pInst->getK(); ++k)
            {
                e1 += q[i][k][t];
                e2 += y[i][k][t];
            }

            std::ostringstream oss;
            oss << "6C_OU" << i << "_" << t;
            constrs.push_back(model.addConstr(e1 >=
                pInst->getUi(i) * e2 - I[i][t], oss.str()));
        }
    }
}


// ok
void init::quantitiesRoutingConstraint(
    GRBModel& model,
    std::vector<GRBConstr>& constrs,
    const std::vector<std::vector<std::vector<GRBVar>>>& y,
    const std::vector<std::vector<std::vector<GRBVar>>>& q,
    const std::shared_ptr<const Instance>& pInst)
{
    DRAW_LOG_F(INFO, "\tinitializing quantities routing constraints");

    for (auto i = 1; i < pInst->getNbVertices(); ++i)
    {
        for (int k = 0; k < pInst->getK(); ++k)
        {
            for (auto t = 0; t < pInst->getT(); ++t)
            {
                std::ostringstream oss;
                oss << "7C" << i << "_" << k << "_" << t;
                constrs.push_back(model.addConstr(q[i][k][t] <=
                                pInst->getUi(i) * y[i][k][t], oss.str()));
            }
        }
    }
}


// ok
void init::capacityVehicleConstraint(
    GRBModel& model,
    std::vector<GRBConstr>& constrs,
    const std::vector<std::vector<std::vector<GRBVar>>>& y,
    const std::vector<std::vector<std::vector<GRBVar>>>& q,
    const std::shared_ptr<const Instance>& pInst)
{
    DRAW_LOG_F(INFO, "\tinitializing capacity vehicles constraints");

    for (int k = 0; k < pInst->getK(); ++k)
    {
        for (auto t = 0; t < pInst->getT(); ++t)
        {
            GRBLinExpr e = 0;
            for (auto i = 1; i < pInst->getNbVertices(); ++i)
            {
                e += q[i][k][t];
            }

            std::ostringstream oss;
            oss << "8C_" << k << "_" << t;
            constrs.push_back(
                model.addConstr(e <= pInst->getCk(k) * y[0][k][t], oss.str()));
        }
    }
}


// ok
void init::degreeConstrs(
    GRBModel& model,
    std::vector<GRBConstr>& constrs,
    const std::vector<std::vector<std::vector<GRBVar>>>& y,
    const std::vector<std::vector<std::vector<std::vector<GRBVar>>>>& x,
    const std::shared_ptr<const Instance>& pInst)
{
    DRAW_LOG_F(INFO, "\tinitializing degree constraints");

    for (auto k = 0; k < pInst->getK(); ++k)
    {
        for (auto t = 0; t < pInst->getT(); ++t)
        {
            for (auto i = 0; i < pInst->getNbVertices(); ++i)
            {
                GRBLinExpr e1 = 0;
                for (auto j = i + 1; j < pInst->getNbVertices(); ++j)
                {
                    e1 += x[i][j][k][t];
                }

                GRBLinExpr e2 = 0;
                for (auto j = 0; j < i; ++j)
                {
                    e2 += x[j][i][k][t];
                }

                std::ostringstream oss;
                oss << "9C_" << i << "_" << k << "_" << t;
                constrs.push_back(
                    model.addConstr(e1 + e2 == 2 * y[i][k][t], oss.str()));
            }
        }
    }
}


void init::subtourEliminationConstrs(
    GRBModel& model,
    std::vector<GRBConstr>& constrs,
    const std::vector<std::vector<std::vector<GRBVar>>>& y,
    const std::vector<std::vector<std::vector<std::vector<GRBVar>>>>& x,
    const std::shared_ptr<const Instance>& pInst)
{
    DRAW_LOG_F(INFO, "\tinitializing subtour elimination constraints");

    const auto nbSets = std::pow(2, pInst->getNbVertices() - 1);
    for (int c = 1; c < nbSets; ++c)
    {
        auto S = findSetS(c, pInst->getNbVertices());
        if (S.size() > 1)
        {
            for (int k = 0; k < pInst->getK(); ++k)
            {
                for (auto t = 0; t < pInst->getT(); ++t)
                {
                    for (auto m : S)
                    {
                        GRBLinExpr lhs = 0;
                        GRBLinExpr rhs = 0;
                        for (auto i : S)
                        {
                            for (auto j : S)
                            {
                                if (i < j)
                                {
                                    lhs += x[i][j][k][t];
                                }
                            }
                            rhs += y[i][k][t];
                        }

                        rhs -= y[m][k][t];

                        std::ostringstream oss;
                        oss << "5CB_" << c << "_" << k << "_" << t << "_" << m;
                        constrs.push_back(model.addConstr(lhs <= rhs,
                                          oss.str()));
                    }
                } 
            }
        }
    }
}