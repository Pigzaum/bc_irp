////////////////////////////////////////////////////////////////////////////////
/*
 * File: init_grb_model.hpp
 * Author: Guilherme O. Chagas
 *
 * @brief Gurobi's model initialization helper functions declartions.
 *
 * (I'm sorry for my bad english xD)
 *
 * Created on October 22, 2020, 10:44 PM.
 * 
 * References:
 */
////////////////////////////////////////////////////////////////////////////////

#ifndef INIT_GRB_MODE_HPP
#define INIT_GRB_MODE_HPP

#include <memory>

#include "gurobi_c++.h"

#include "instance.hpp"

namespace init
{

/**
 * @brief:.
 * @param:.
 * @param:.
 * @param:.
*/
void inventoryLevelVariables(GRBModel& model,
                             std::vector<std::vector<GRBVar>>& I,
                             const std::shared_ptr<const Instance>& pInst);

/**
 * @brief:.
 * @param:.
 * @param:.
 * @param:.
*/
void quantityVariables(GRBModel& model,
                       std::vector<std::vector<GRBVar>>& q,
                       const std::shared_ptr<const Instance>& pInst);


/**
 * @brief:.
 * @param:.
 * @param:.
 * @param:.
*/
void visitationVariables(GRBModel& model,
                         std::vector<std::vector<GRBVar>>& y,
                         const std::shared_ptr<const Instance>& pInst);

/**
 * @brief:.
 * @param:.
 * @param:.
 * @param:.
*/
void routingVariables(GRBModel& model,
                      std::vector<std::vector<std::vector<GRBVar>>>& x,
                      const std::shared_ptr<const Instance>& pInst);

/**
 * @brief:.
 * @param:.
 * @param:.
 * @param:.
 * @param:.
 * @param:.
*/
void inventoryDefDepotConstrs(GRBModel& model,
                              std::vector<GRBConstr>& constrs,
                              const std::vector<std::vector<GRBVar>>& I,
                              const std::vector<std::vector<GRBVar>>& q,
                              const std::shared_ptr<const Instance>& pInst);

/**
 * @brief:.
 * @param:.
 * @param:.
 * @param:.
 * @param:.
 * @param:.
*/
void stockOutDepotConstrs(GRBModel& model,
                          std::vector<GRBConstr>& constrs,
                          const std::vector<std::vector<GRBVar>>& I,
                          const std::vector<std::vector<GRBVar>>& q,
                          const std::shared_ptr<const Instance>& pInst);

/**
 * @brief:.
 * @param:.
 * @param:.
 * @param:.
 * @param:.
 * @param:.
*/
void inventoryDefCustomersConstrs(GRBModel& model,
                                  std::vector<GRBConstr>& constrs,
                                  const std::vector<std::vector<GRBVar>>& I,
                                  const std::vector<std::vector<GRBVar>>& q,
                                  const std::shared_ptr<const Instance>& pInst);

/**
 * @brief:.
 * @param:.
 * @param:.
 * @param:.
 * @param:.
 * @param:.
*/
void capacityConstrs(GRBModel& model,
                     std::vector<GRBConstr>& constrs,
                     const std::vector<std::vector<GRBVar>>& q,
                     const std::shared_ptr<const Instance>& pInst);

/**
 * @brief:.
 * @param:.
 * @param:.
 * @param:.
 * @param:.
 * @param:.
*/
void inventoryLevelConstrs(GRBModel& model,
                           std::vector<GRBConstr>& constrs,
                           const std::vector<std::vector<GRBVar>>& I,
                           const std::vector<std::vector<GRBVar>>& q,
                           const std::shared_ptr<const Instance>& pInst);

/**
 * @brief:.
 * @param:.
 * @param:.
 * @param:.
 * @param:.
 * @param:.
*/
void mlQuantityCapacityConstrs(GRBModel& model,
                               std::vector<GRBConstr>& constrs,
                               const std::vector<std::vector<GRBVar>>& I,
                               const std::vector<std::vector<GRBVar>>& q,
                               const std::shared_ptr<const Instance>& pInst);

/**
 * @brief:.
 * @param:.
 * @param:.
 * @param:.
 * @param:.
 * @param:.
*/
void ouQuantityCapacityConstrs(GRBModel& model,
                               std::vector<GRBConstr>& constrs,
                               const std::vector<std::vector<GRBVar>>& I,
                               const std::vector<std::vector<GRBVar>>& y,
                               const std::vector<std::vector<GRBVar>>& q,
                               const std::shared_ptr<const Instance>& pInst);

/**
 * @brief:.
 * @param:.
 * @param:.
 * @param:.
 * @param:.
 * @param:.
*/
void quantitiesRoutingConstraint(GRBModel& model,
                                 std::vector<GRBConstr>& constrs,
                                 const std::vector<std::vector<GRBVar>>& y,
                                 const std::vector<std::vector<GRBVar>>& q,
                                 const std::shared_ptr<const Instance>& pInst);

/**
 * @brief:.
 * @param:.
 * @param:.
 * @param:.
 * @param:.
 * @param:.
*/
void capacityVehicleConstraint(GRBModel& model,
                               std::vector<GRBConstr>& constrs,
                               const std::vector<std::vector<GRBVar>>& y,
                               const std::vector<std::vector<GRBVar>>& q,
                               const std::shared_ptr<const Instance>& pInst);

/**
 * @brief:.
 * @param:.
 * @param:.
 * @param:.
 * @param:.
 * @param:.
*/
void degreeConstrs(GRBModel& model,
                   std::vector<GRBConstr>& constrs,
                   const std::vector<std::vector<GRBVar>>& y,
                   const std::vector<std::vector<std::vector<GRBVar>>>& x,
                   const std::shared_ptr<const Instance>& pInst);


/**
 * @brief:.
 * @param:.
 * @param:.
 * @param:.
 * @param:.
 * @param:.
*/
void subtourEliminationConstrs(
    GRBModel& model,
    std::vector<GRBConstr>& constrs,
    const std::vector<std::vector<GRBVar>>& y,
    const std::vector<std::vector<std::vector<GRBVar>>>& x,
    const std::shared_ptr<const Instance>& pInst);

} // init namespace

#endif // INIT_IRP_MODE_HPP