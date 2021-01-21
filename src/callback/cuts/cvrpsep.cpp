////////////////////////////////////////////////////////////////////////////////
/*
 * File: cvrpsep.cpp
 * Author: Guilherme O. Chagas
 *
 * @brief Lazy and cut CVRPSEP functions.
 *
 * (I'm sorry for my bad english xD)
 *
 * Created on January 20, 2020, 06:51 PM
 * 
 * References:
 */
////////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "../../../include/ext/cvrpsep/capsep.h"
#include "../../../include/ext/cvrpsep/cnstrmgr.h"
#include "../../../include/ext/loguru/loguru.hpp"

#include "../../../include/callback/callback_sec.hpp"
#include "../../../include/utils/constants.hpp"

////////////////////////////// Helper functions  ///////////////////////////////

namespace
{

static const int cDim = 100; // ?
static const int cMaxNbCapCuts = 8; // ?

/**
 * @brief
*/
int checkForDepot(const int i, const int n)
{
    return i == n ? 0 : i;
}

/*
* @brief
*/
std::vector<int> getDemand(
    const std::vector<std::vector<std::vector<double>>>& qVal,
    const int k,
    const int t)
{
    std::vector<int> demand(qVal.size(), 0);

    for (std::size_t i = 1; i < demand.size(); ++i)
    {
        demand[i] = static_cast<int>(std::ceil(qVal[i][k][t]));
    }

    return demand;
}

/**
 * @brief
*/
std::tuple<int, std::vector<int>, std::vector<int>, std::vector<double>>
    getEdges(
        const std::vector<std::vector<std::vector<std::vector<double>>>>& xVal,
        const int n,
        const int k,
        const int t)
{
    int nbEdges = 0;
    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            if (xVal[i][j][k][t] > utils::GRB_EPSILON)
            {
                ++nbEdges;
            }
        }
    }
    
    std::vector<int> edgeTail;
    std::vector<int> edgeHead;
    std::vector<double> edgeX;

    if (nbEdges > 0)
    {
        edgeTail.reserve(nbEdges + 1);
        edgeHead.reserve(nbEdges + 1);
        edgeX.reserve(nbEdges + 1);
        /* skip first element */
        edgeTail.push_back(0); 
        edgeHead.push_back(0);
        edgeX.push_back(0);

        for (int i = 0; i < n; ++i)
        {
            for (int j = i + 1; j < n; ++j)
            {
                if (xVal[i][j][k][t] > utils::GRB_EPSILON)
                {
                    edgeTail.push_back(i == 0 ? n : i);
                    edgeHead.push_back(j);
                    edgeX.push_back(xVal[i][j][k][t]);
                }
            }
        }
    }

    return std::make_tuple(nbEdges, edgeTail, edgeHead, edgeX);
}

}

////////////////////////////////////////////////////////////////////////////////


int CallbackSEC::addCVRPSEPCAP(const constrsType cstType)
{
    int nbAdded = 0; // # of lazy/cuts added in this call

    /* get solution: delivery quantities */
    auto qVal = getqVarsValues(cstType);
    /* get solution: routing */
    auto xVal = getxVarsValues(cstType);
    /* get solution: visitation */
    auto yVal = getyVarsValues(cstType);

    const int n = mpInst->getNbVertices();
    for (int k = 0; k < mpInst->getK(); ++k)
    {
        for (int t = 0; t < mpInst->getT(); ++t)
        {
            auto [nbEdges, edgeTail, edgeHead, edgeX] = getEdges(xVal, n, k, t);

            if (nbEdges == 0) continue;

            /* Parameters of the CVRPSEP */
            char integerAndFeasible;
            double maxViolation = 0;
            auto demand = getDemand(qVal, k, t);
            CnstrMgrPointer cutsCMP, myOldCutsCMP;
            CMGR_CreateCMgr(&cutsCMP, cDim);
            CMGR_CreateCMgr(&myOldCutsCMP, cDim);

            CAPSEP_SeparateCapCuts(n - 1,
                                   demand.data(),
                                   mpInst->getCk(k),
                                   nbEdges,
                                   edgeTail.data(),
                                   edgeHead.data(),
                                   edgeX.data(),
                                   myOldCutsCMP,
                                   cMaxNbCapCuts,
                                   utils::GRB_EPSILON,
                                   &integerAndFeasible,
                                   &maxViolation,
                                   cutsCMP);

            for (int c = 0; c < cutsCMP->Size; ++c)
            {
                if (cutsCMP->CPL[c]->CType == CMGR_CT_CAP)
                {
                    std::vector<int> list;
                    list.reserve(cutsCMP->CPL[c]->IntListSize);
                    for (int j = 1; j <= cutsCMP->CPL[c]->IntListSize; ++j)
                    {
                        list.push_back(
                            checkForDepot(cutsCMP->CPL[c]->IntList[j], n));
                    }

                    GRBLinExpr xExpr = 0, yExpr = 0;
                    double xSum = 0, ySum = 0;
                    for (std::size_t i = 0; i < list.size(); ++i)
                    {
                        for (std::size_t j = 0; j < list.size(); ++j)
                        {
                            if (list[i] < list[j])
                            {
                                xExpr += m_x[list[i]][list[j]][k][t];
                                xSum += xVal[list[i]][list[j]][k][t];
                            }
                        }
                        yExpr += m_y[list[i]][k][t];
                        ySum += yVal[list[i]][k][t];
                    }

                    for (std::size_t i = 0; i < list.size(); ++i)
                    {
                        if (cstType == constrsType::lazy)
                        {
                            addLazy(xExpr <= yExpr - m_y[list[i]][k][t]);
                            ++nbAdded;
                        }
                        else if (xSum - ySum + yVal[list[i]][k][t] >
                                 utils::GRB_EPSILON)
                        {
                            addCut(xExpr <= yExpr - m_y[list[i]][k][t]);
                            ++nbAdded;
                        }
                    }
                }
            }

            CMGR_FreeMemCMgr(&cutsCMP);
            CMGR_FreeMemCMgr(&myOldCutsCMP);
        }
    }

    return nbAdded;
}