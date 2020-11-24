////////////////////////////////////////////////////////////////////////////////
/*
 * File: instance.hpp
 * Author: Guilherme O. Chagas
 *
 * @brief IRP instance class declaration.
 *
 * (I'm sorry for my bad english xD)
 *
 * Created on October 22, 2020, 06:44 PM.
 * 
 * References:
 * https://bit.ly/3dQu0Kc
 */
////////////////////////////////////////////////////////////////////////////////

#ifndef INSTANCE_HPP
#define INSTANCE_HPP

#include <string>
#include <vector>

class ConfigParameters;

class Instance
{
public:

    /**
     * @brief Default constructor, copy constructor, move constructor,
     * destructor, copy assingment operator and move assingment operator.
    */
    Instance() = default;
    Instance(const Instance& other) = default;
    Instance(Instance&& other) = default;
    ~Instance() = default;
    Instance& operator=(const Instance& other) = default;
    Instance& operator=(Instance&& other) = default;

    /**
     * @brief Constructs from a instance file.
     * @param : const std::string&: instance file path.
    */
    Instance(const std::string& filePath, const int K);

    double getC() const;

    double getCk(const int k) const;

    double get_cij(const int i, const int j) const;

    double get_hi(const int i) const;

    double getIi0(const int i) const;

    int getK() const; // number of vehicles

    double getLi(const int i) const;

    std::string getName() const;

    int getNbVertices() const;

    double get_rit(const int i, const int t) const;

    int getT() const;

    double getUi(const int i) const;

    void setK(const int K);

    /**
     * @brief: Print instance on console.
    */
    void show() const;

private:

    // instance full path
    std::string mPath;

    // number of vertices (depot and customers)
    int mNbVertices; 

    // number of vehicles
    int mK;

    // number of discrete time instants of the planning time horizon
    int mT;

    // transportation capacity
    double mC;

    // transportation capacity of vehicle k
    std::vector<int> mCk;

    // unit inventory cost
    std::vector<double> m_hi;

    // starting level of the inventory
    std::vector<double> mIi0;

    // minimum level of the inventory at the retailer i
    std::vector<double> mLi;

    // maximum level of the inventory at the retailer i
    std::vector<double> mUi;

    // quantity absorved by the retailer i at each dicrete time instante
    std::vector<std::vector<double>> m_rit;

    // (x, y) coordinates of the supply and the retailers
    std::vector<std::pair<double, double>> mCoord;

    // distance matrix
    std::vector<std::vector<int>> m_cij;

    /**
     * @brief Initializes the instance object from the instance file path.
     * @param:.
    */
    void init(const std::string& filePath);
};

#endif // INSTANCE_HPP