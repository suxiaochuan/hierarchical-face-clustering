//
// Created by Haowen Xu on 2018/6/6.
//

#ifndef HIERARCHICALFACECLUSTERING_CLUSTERINGALGORITHM_H
#define HIERARCHICALFACECLUSTERING_CLUSTERINGALGORITHM_H

#include <unordered_set>
#include <Eigen/Dense>

struct DualNodeType {
  std::unordered_set<int> faces;
  Eigen::Matrix3d A;
  Eigen::Vector3d b;
  double c;

  DualNodeType() = default;
  DualNodeType(int f, const Eigen::Matrix3d &A, const Eigen::Vector3d &b, double c);

  DualNodeType operator+(const DualNodeType &rhs) const;
};

class ClusteringAlgorithm {
 public:
  // constructors
  ClusteringAlgorithm() = default;

  std::vector<std::vector<int>> Cluster(const std::string &input_filename, double th);

 private:

  void GetDualEdges(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, std::vector<std::pair<int, int>> &dual_edges);
  double EvaluateCost(const DualNodeType &dn1, const DualNodeType &dn2);
};

#endif //HIERARCHICALFACECLUSTERING_CLUSTERINGALGORITHM_H
