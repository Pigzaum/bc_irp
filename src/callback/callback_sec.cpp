////////////////////////////////////////////////////////////////////////////////
/*
 * File: callback_sec.cpp
 * Author: Guilherme O. Chagas
 *
 * @brief Callback class definition for lazy/cut subtour separation
 * constraints.
 *
 * @acknowledgment Special thanks to Ph.D. Cleder Marcos Schenekemberg.
 *
 * (I'm sorry for my bad english xD)
 *
 * Created on November 16, 2020, 00:15 AM
 * 
 * References:
 */
////////////////////////////////////////////////////////////////////////////////

#include "../../include/ext/loguru/loguru.hpp"

#include "../../include/callback/callback_sec.hpp"
#include "../../include/ext/cvrpsep/capsep.h"
#include "../../include/ext/cvrpsep/cnstrmgr.h"

////////////////////////////////////////////////////////////////////////////////

namespace
{

static const int C_EPS = 1e-5;

} // anonymous namespace

////////////////////////////////////////////////////////////////////////////////

CallbackSEC::CallbackSEC(
    const std::vector<std::vector<std::vector<std::vector<GRBVar>>>>& x,
    const std::vector<std::vector<std::vector<GRBVar>>>& y,
    const std::shared_ptr<const Instance>& p_inst) :
        m_x(x),
        m_y(y),
        mp_inst(p_inst)
{}


void CallbackSEC::callback()
{
    try
    {
        switch (where)
        {
        case GRB_CB_MIPSOL:
        {
            addLazyCVRPSEP();
            break;
        }
        case GRB_CB_MIPNODE:
        {
            if (getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL)
            {
                addCutCVRPSEP();
            }
            break;
        }
        default:
            break;
        }
    }
    catch (GRBException& e)
    {
        RAW_LOG_F(ERROR, "callback() exp: %s", e.getMessage().c_str());
    }
    catch (...)
    {
        RAW_LOG_F(ERROR, "callback(): Unknown Exception");
    }
}


void CallbackSEC::addLazyCVRPSEP()
{
    std::vector<std::vector<std::vector<std::vector<double>>>>
        valueX(mp_inst->getNbVertices(),
             std::vector<std::vector<std::vector<double>>>(
                mp_inst->getNbVertices(),
                std::vector<std::vector<double>>(
                    mp_inst->getK(),
                    std::vector<double>(mp_inst->getT(), 0))));

    /* get solution: routing */
    for (int i = 0; i < mp_inst->getNbVertices(); ++i)
    {
        for (int j = i + 1; j < mp_inst->getNbVertices(); ++j)
        {
            for (int k = 0; k < mp_inst->getK(); ++k)
            {
                for (int t = 0; t < mp_inst->getT(); ++t)
                {
                    valueX[i][j][k][t] = getSolution(m_x[i][j][k][t]);
                }
            }
        }
    }

    for (int k = 0; k < mp_inst->getK(); ++k)
    {
        for (int t = 0; t < mp_inst->getT(); ++t)
        {
            int nbEdges = 0;
            for (int i = 0; i < mp_inst->getNbVertices(); ++i)
            {
                for (int j = i + 1; j < mp_inst->getNbVertices(); ++j)
                {
                    if (valueX[i][j][k][t] > C_EPS)
                    {
                        ++nbEdges;
                    }
                }
            }

            if (nbEdges > 0)
            {
                /* Parameters of the CVRPSEP */
                const int maxNbCuts = 8;
                char integerAndFeasible;
                double maxViolation = 0;

                std::vector<int> demand(mp_inst->getNbVertices());

                for (int i = 1; i < mp_inst->getNbVertices(); ++i)
                {
                    demand[i] = mp_inst->get_rit(i, t);
                }

                std::vector<int> edgeTail(nbEdges + 1, 0);
                std::vector<int> edgeHead(nbEdges + 1, 0);
                std::vector<double> edgeX(nbEdges + 1, 0);

                int aux = 1;
                for (int i = 0; i < mp_inst->getNbVertices(); ++i)
                {
                    for (int j = i + 1; j < mp_inst->getNbVertices(); ++j)
                    {
                        if (valueX[i][j][k][t] > C_EPS)
                        {
                            if (i == 0)
                            {
                                edgeTail[aux] = mp_inst->getNbVertices();
                            }
                            else
                            {
                                edgeTail[aux] = i;
                            }
                            edgeHead[aux] = j;
                            edgeX[aux] = valueX[i][j][k][t];
                            ++aux;
                        }
                    }
                }

                const int dim = 100;
                CnstrMgrPointer myCutsCMP, myOldCutsCMP;
                CMGR_CreateCMgr(&myCutsCMP, dim);
                CMGR_CreateCMgr(&myOldCutsCMP, dim);

                CAPSEP_SeparateCapCuts(mp_inst->getNbVertices() - 1,
                                    demand.data(),
                                    mp_inst->getCk(k),
                                    nbEdges,
                                    edgeTail.data(),
                                    edgeHead.data(),
                                    edgeX.data(),
                                    myOldCutsCMP,
                                    maxNbCuts,
                                    C_EPS,
                                    &integerAndFeasible,
                                    &maxViolation,
                                    myCutsCMP);

                /* capacity */
                std::vector<int> list(mp_inst->getNbVertices(), 0);

                for (int cut = 0; cut < myCutsCMP->Size; ++cut)
                {
                    if (myCutsCMP->CPL[cut]->CType == CMGR_CT_CAP)
                    {
                        /* capacity cut */
                        int listSize = 0;
                        for (int j = 1;
                             j <= myCutsCMP->CPL[cut]->IntListSize; ++j)
                        {
                            list[++listSize] = myCutsCMP->CPL[cut]->IntList[j];
                        }

                        GRBLinExpr e1 = 0;
                        for (int i = 1; i <= listSize; ++i)
                        {
                            for (int j = 1; j <= listSize; ++j)
                            {
                                if (list[i] < list[j])
                                {
                                    e1 += m_x[list[i]][list[j]][k][t];
                                }
                            }
                        }

                        GRBLinExpr e2 = 0;
                        for (int i = 1; i <= listSize; ++i)
                        {
                            e2 += m_y[list[i]][k][t];
                        }

                        for (int i = 1; i <= listSize; ++i)
                        {
                            addLazy(e1 <= e2 - m_y[list[i]][k][t]);
                        }
                    }
                }

                CMGR_FreeMemCMgr(&myCutsCMP);
                CMGR_FreeMemCMgr(&myOldCutsCMP);
            }
        }
    }
}


void CallbackSEC::addCutCVRPSEP()
{
    std::vector<std::vector<std::vector<double>>> valueY(
        mp_inst->getNbVertices(),
        std::vector<std::vector<double>>(
            mp_inst->getK(),
            std::vector<double>(mp_inst->getT(), 0)));

    /* get solution: routing */
    for (int i = 0; i < mp_inst->getNbVertices(); ++i)
    {
        for (int k = 0; k < mp_inst->getK(); ++k)
        {
            for (int t = 0; t < mp_inst->getT(); ++t)
            {
                valueY[i][k][t] = getNodeRel(m_y[i][k][t]);
            }
        }
    }

    std::vector<std::vector<std::vector<std::vector<double>>>>
        valueX(mp_inst->getNbVertices(),
             std::vector<std::vector<std::vector<double>>>(
                mp_inst->getNbVertices(),
                std::vector<std::vector<double>>(
                    mp_inst->getK(),
                    std::vector<double>(mp_inst->getT(), 0))));

    /* get solution: routing */
    for (int i = 0; i < mp_inst->getNbVertices(); ++i)
    {
        for (int j = i + 1; j < mp_inst->getNbVertices(); ++j)
        {
            for (int k = 0; k < mp_inst->getK(); ++k)
            {
                for (int t = 0; t < mp_inst->getT(); ++t)
                {
                    valueX[i][j][k][t] = getNodeRel(m_x[i][j][k][t]);
                }
            }
        }
    }

    for (int k = 0; k < mp_inst->getK(); ++k)
    {
        for (int t = 0; t < mp_inst->getT(); ++t)
        {
            int nbEdges = 0;
            for (int i = 0; i < mp_inst->getNbVertices(); ++i)
            {
                for (int j = i + 1; j < mp_inst->getNbVertices(); ++j)
                {
                    if (valueX[i][j][k][t] > C_EPS)
                    {
                        ++nbEdges;
                    }
                }
            }

            if (nbEdges > 0)
            {
                /* Parameters of the CVRPSEP */
                std::vector<int> demand(mp_inst->getNbVertices(), 0);
                std::vector<int> edgeTail(nbEdges + 1);
                std::vector<int> edgeHead(nbEdges + 1);
                std::vector<double> edgeX(nbEdges + 1);

                char integerAndFeasible;
                double maxViolation = 0;

                const int maxNbCapCuts = 8;

                for (int i = 1; i < mp_inst->getNbVertices(); ++i)
                {
                    demand[i] = mp_inst->get_rit(i, t);
                }

                int aux = 1;
                for (int i = 0; i < mp_inst->getNbVertices(); ++i)
                {
                    for (int j = i + 1; j < mp_inst->getNbVertices(); ++j)
                    {
                        if (valueX[i][j][k][t] > C_EPS)
                        {
                            if (i == 0)
                            {
                                edgeTail[aux] = mp_inst->getNbVertices();
                            }
                            else
                            {
                                edgeTail[aux] = i;
                            }
                            edgeHead[aux] = j;
                            edgeX[aux] = valueX[i][j][k][t];
                            ++aux;
                        }
                    }
                }

                const int dim = 100;
                CnstrMgrPointer myCutsCMP, myOldCutsCMP;
                CMGR_CreateCMgr(&myCutsCMP, dim);
                CMGR_CreateCMgr(&myOldCutsCMP, dim);

                CAPSEP_SeparateCapCuts(mp_inst->getNbVertices() - 1,
                                    demand.data(),
                                    mp_inst->getCk(k),
                                    nbEdges,
                                    edgeTail.data(),
                                    edgeHead.data(),
                                    edgeX.data(),
                                    myOldCutsCMP,
                                    maxNbCapCuts,
                                    C_EPS,
                                    &integerAndFeasible,
                                    &maxViolation,
                                    myCutsCMP);

                /* capacity */
                std::vector<int> list(mp_inst->getNbVertices(), 0);

                for (int cut = 0; cut < myCutsCMP->Size; ++cut)
                {
                    if (myCutsCMP->CPL[cut]->CType == CMGR_CT_CAP)
                    {
                        /* capacity cut */
                        int listSize = 0;
                        for (int j = 1;
                             j <= myCutsCMP->CPL[cut]->IntListSize; ++j)
                        {
                            list[++listSize] = myCutsCMP->CPL[cut]->IntList[j];
                        }

                        //double rhs = myCutsCMP->CPL[cut]->RHS;
                        GRBLinExpr e1 = 0;
                        double sumX = 0;
                        for (int i = 1; i <= listSize; ++i)
                        {
                            for (int j = 1; j <= listSize; ++j)
                            {
                                if (list[i] < list[j])
                                {
                                    e1 += m_x[list[i]][list[j]][k][t];
                                    sumX += valueX[list[i]][list[j]][k][t];
                                }
                            }
                        }

                        GRBLinExpr e2 = 0;
                        double sumY = 0;
                        for (int i = 1; i <= listSize; ++i)
                        {
                            e2 += m_y[list[i]][k][t];
                            sumY += valueY[list[i]][k][t];
                        }

                        for (int i = 1; i <= listSize; ++i)
                        {
                            if (sumX >= (sumY - valueY[list[i]][k][t]))
                            {
                                addCut(e1 <= e2 - m_y[list[i]][k][t]);
                            }
                        }
                    }
                }

                CMGR_FreeMemCMgr(&myCutsCMP);
                CMGR_FreeMemCMgr(&myOldCutsCMP);
            }
        }
    }
}