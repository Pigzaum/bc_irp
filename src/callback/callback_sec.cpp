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
    const std::vector<std::vector<std::vector<GRBVar>>>& q,
    const std::vector<std::vector<std::vector<std::vector<GRBVar>>>>& x,
    const std::vector<std::vector<std::vector<GRBVar>>>& y,
    const std::shared_ptr<const Instance>& p_inst) :
        m_q(q),
        m_x(x),
        m_y(y),
        mpInst(p_inst)
{}


void CallbackSEC::callback()
{
    try
    {
        if (where == GRB_CB_MIPSOL)
        {
            addCVRPSEPCAP(constrsType::lazy);
        }
        else if (where == GRB_CB_MIPNODE &&
                 getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL)
        {
            addCVRPSEPCAP(constrsType::cut);
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

/* -------------------------------------------------------------------------- */

std::vector<std::vector<std::vector<double>>> CallbackSEC::getqVarsValues(
    const constrsType cstType)
{
    std::vector<std::vector<std::vector<double>>> qVal(
        mpInst->getNbVertices(), std::vector<std::vector<double>>(
            mpInst->getK(), std::vector<double>(mpInst->getT(), 0)
        )
    );

    for (int i = 1; i < mpInst->getNbVertices(); ++i)
    {
        for (int k = 0; k < mpInst->getK(); ++k)
        {
            for (int t = 0; t < mpInst->getT(); ++t)
            {
                DCHECK_F(i < static_cast<int>(m_q.size()));
                DCHECK_F(k < static_cast<int>(m_q[i].size()));
                DCHECK_F(t < static_cast<int>(m_q[i][k].size()));
                if (cstType == constrsType::lazy)
                {
                    qVal[i][k][t] = getSolution(m_q[i][k][t]);
                }
                else
                {
                    qVal[i][k][t] = getNodeRel(m_q[i][k][t]);
                }
            }
        }
    }

    return qVal;
}


std::vector<std::vector<std::vector<std::vector<double>>>>
    CallbackSEC::getxVarsValues(const constrsType cstType)
{
    std::vector<std::vector<std::vector<std::vector<double>>>> xVal(
        mpInst->getNbVertices(), std::vector<std::vector<std::vector<double>>>(
            mpInst->getNbVertices(), std::vector<std::vector<double>>(
                mpInst->getK(), std::vector<double>(mpInst->getT(), 0)
            )
        )
    );

    for (int i = 0; i < mpInst->getNbVertices(); ++i)
    {
        for (int j = i + 1; j < mpInst->getNbVertices(); ++j)
        {
            for (int k = 0; k < mpInst->getK(); ++k)
            {
                for (int t = 0; t < mpInst->getT(); ++t)
                {
                    DCHECK_F(i < static_cast<int>(m_x.size()));
                    DCHECK_F(j < static_cast<int>(m_x[i].size()));
                    DCHECK_F(k < static_cast<int>(m_x[i][j].size()));
                    DCHECK_F(t < static_cast<int>(m_x[i][j][k].size()));
                    if (cstType == constrsType::lazy)
                    {
                        xVal[i][j][k][t] = getSolution(m_x[i][j][k][t]);
                    }
                    else
                    {
                        xVal[i][j][k][t] = getNodeRel(m_x[i][j][k][t]);
                    }
                }
            }
        }
    }

    return xVal;
}


std::vector<std::vector<std::vector<double>>> CallbackSEC::getyVarsValues(
        const constrsType cstType)
{
    std::vector<std::vector<std::vector<double>>> yVal(
        mpInst->getNbVertices(), std::vector<std::vector<double>>(
            mpInst->getK(), std::vector<double>(mpInst->getT(), 0)
        )
    );

    for (int i = 0; i < mpInst->getNbVertices(); ++i)
    {
        for (int k = 0; k < mpInst->getK(); ++k)
        {
            for (int t = 0; t < mpInst->getT(); ++t)
            {
                DCHECK_F(i < static_cast<int>(m_y.size()));
                DCHECK_F(k < static_cast<int>(m_y[i].size()));
                DCHECK_F(t < static_cast<int>(m_y[i][k].size()));
                if (cstType == constrsType::lazy)
                {
                    yVal[i][k][t] = getSolution(m_y[i][k][t]);
                }
                else
                {
                    yVal[i][k][t] = getNodeRel(m_y[i][k][t]);
                }
            }
        }
    }

    return yVal;
}