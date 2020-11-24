////////////////////////////////////////////////////////////////////////////////
/*
 * File: instance.cpp
 * Author: Guilherme O. Chagas
 *
 * @brief IRP instance class declaration.
 *
 * (I'm sorry for my bad english xD)
 *
 * Created on October 22, 2020, 07:02 PM.
 * 
 * References:
 * https://bit.ly/3dQu0Kc
 */
////////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>

#include "../include/ext/loguru/loguru.hpp"

#include "../include/config_parameters.hpp"
#include "../include/instance.hpp"

/////////////////////////////// Helper methods  ////////////////////////////////

namespace
{

/**
 * @brief:.
 * @param:.
 * @return:.
*/
std::vector<std::vector<int>> compute_distance_mtx(
    const std::vector<std::pair<double, double>>& coord)
{
    std::vector<std::vector<int>> cij(coord.size(),
                                      std::vector<int>(coord.size(), 0));

    for (std::size_t i = 0; i < coord.size(); ++i)
    {
        for (std::size_t j = i + 1; j < coord.size(); ++j)
        {
            auto tmp = std::sqrt(std::pow(coord[i].first - coord[j].first, 2) +
                                 std::pow(coord[i].second - coord[j].second,2));
            cij[i][j] = static_cast<int>(std::round(tmp));
            cij[j][i] = cij[i][j];
        }
    }

    return cij;
}


template <typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& v)
{
    if (!v.empty())
    {
        out << '[';
        std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
        out << "\b\b]";
    }
    return out;
}


} // anonymous namespace

////////////////////////////////////////////////////////////////////////////////

Instance::Instance(const std::string& filePath, const int K) :
    mPath(filePath),
    mK(K)
{
    CHECK_F(std::filesystem::exists(filePath));
    init(filePath);
}


double Instance::getC() const
{
    return mC;
}


double Instance::getCk(const int k) const
{
    DCHECK_F(k < static_cast<int>(mCk.size()));
    return mCk[k];
}


double Instance::get_cij(const int i, const int j) const
{
    DCHECK_F(i < static_cast<int>(m_cij.size()));
    DCHECK_F(j < static_cast<int>(m_cij[i].size()));
    return m_cij[i][j];
}


double Instance::get_hi(const int i) const
{
    DCHECK_F(i < static_cast<int>(m_hi.size()));
    return m_hi[i];
}


double Instance::getIi0(const int i) const
{
    DCHECK_F(i < static_cast<int>(mIi0.size()));
    return mIi0[i];
}


double Instance::getLi(const int i) const
{
    DCHECK_F(i < static_cast<int>(mLi.size()));
    return mLi[i];
}


int Instance::getK() const
{
    return mK;
}


std::string Instance::getName() const
{
    return std::filesystem::path(mPath).stem().string();
}


int Instance::getNbVertices() const
{
    return mNbVertices;
}


double Instance::get_rit(const int i, const int t) const
{
    DCHECK_F(i < static_cast<int>(m_rit.size()));
    DCHECK_F(t < static_cast<int>(m_rit[i].size()));
    return m_rit[i][t];
}


int Instance::getT() const
{
    return mT;
}


double Instance::getUi(const int i) const
{
    DCHECK_F(i < static_cast<int>(mUi.size()));
    return mUi[i];
}


void Instance::setK(const int K)
{
    mK = K;
}


void Instance::show() const
{
    RAW_LOG_F(INFO, std::string(80, '-').c_str());

    std::ostringstream oss;

    oss << "Instance details\n";
    oss << "n\t" << mNbVertices << "\n";
    oss << "H\t" << mT << "\n";
    oss << "C\t" << mC;

    RAW_LOG_F(INFO, oss.str().c_str());

    oss.clear();
    oss.str("");
    oss << "Li:\n";
    oss << mLi << "\n";
    oss << "Ui:\n";
    oss << mUi << "\n";

    oss << "Distance matrix:\n";
    for (auto &r : m_cij)
    {
        oss << "|";
        std::string sep = "";
        for (auto e : r)
        {
            oss << sep << e;
            sep = "\t";
        }
        oss << sep << "|\n";
    }

    RAW_LOG_F(INFO, oss.str().c_str());
    RAW_LOG_F(INFO, std::string(80, '-').c_str());
}


/////////////////////////////// private methods ////////////////////////////////

void Instance::init(const std::string& filePath)
{
    std::ifstream file(filePath);

    file >> mNbVertices;
    file >> mT;
    file >> mC;

    m_hi.reserve(mNbVertices);
    mIi0.reserve(mNbVertices);
    m_rit.reserve(mNbVertices);
    mLi.reserve(mNbVertices);
    mLi.push_back(0);
    mUi.reserve(mNbVertices);
    mUi.push_back(0);
    mCoord.reserve(mNbVertices);

    mCk.reserve(mK);
    for (int i = 0; i < mK; ++i)
    {
        mCk.push_back(std::round(mC / mK));
    }

    /* first line: depot */
    {
        double x, y;
        file >> x >> x >> y; // skip index
        mCoord.push_back(std::make_pair(x, y));
    }

    double tmp = 0;
    file >> tmp;
    mIi0.push_back(tmp);

    file >> tmp;
    m_rit = std::vector<std::vector<double>>(mNbVertices,
                                             std::vector<double>(mT + 1, 0));
    for (auto t = 0; t < mT; ++t)
    {
        DCHECK_F(t < static_cast<int>(m_rit[0].size()));
        m_rit[0][t] = tmp;
    }

    file >> tmp;
    m_hi.push_back(tmp);

    /* retailers */
    std::string line;
    int count = 1;
    while (count < mNbVertices)
    {
        int index = 0;
        file >> index;

        {
            double x, y;
            file >> x >> y;
            mCoord.push_back(std::make_pair(x, y));
        }

        file >> tmp;
        mIi0.push_back(tmp);
        file >> tmp;
        mUi.push_back(tmp);
        file >> tmp;
        mLi.push_back(tmp);

        file >> tmp;
        for (auto t = 0; t < mT; ++t)
        {
            DCHECK_F(index - 1 < static_cast<int>(m_rit.size()));
            DCHECK_F(t < static_cast<int>(m_rit[index - 1].size()));
            m_rit[index - 1][t] = tmp;
        }

        file >> tmp;
        m_hi.push_back(tmp);
        ++count;
    }

    m_cij = compute_distance_mtx(mCoord);
}