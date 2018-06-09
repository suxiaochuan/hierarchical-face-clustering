//
// Created by Haowen Xu on 2018/6/6.
//

#include "ClusteringAlgorithm.h"
#include <queue>
#include <algorithm>
#include <igl/readOBJ.h>
#include <igl/edge_topology.h>


DualNodeType::DualNodeType(int f, const Eigen::Matrix3d &A, const Eigen::Vector3d &b, double c)  : A(A), b(b), c(c) {
  faces = std::unordered_set<int>();
  faces.insert(f);
}

DualNodeType DualNodeType::operator+(const DualNodeType &rhs) const {
  DualNodeType res;
  res.faces = faces;
  for (int face : rhs.faces)
    res.faces.insert(face);
  res.A = A + rhs.A;
  res.b = b + rhs.b;
  res.c = c + rhs.c;
  return res;
}

std::vector<std::vector<int>> ClusteringAlgorithm::Cluster(const std::string &input_filename, double th) {
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::readOBJ(input_filename, V, F);

  std::vector<DualNodeType> dual_nodes;
  std::unordered_set<int> valid;
  for (int i = 0; i < F.rows(); ++i) {
    int i0 = F(i, 0), i1 = F(i, 1), i2 = F(i, 2);
    DualNodeType dn0(i, V.row(i0).transpose().eval() * V.row(i0), V.row(i0).transpose(), 1);
    DualNodeType dn1(i, V.row(i1).transpose().eval() * V.row(i1), V.row(i1).transpose(), 1);
    DualNodeType dn2(i, V.row(i2).transpose().eval() * V.row(i2), V.row(i2).transpose(), 1);
    dual_nodes.emplace_back(dn0 + dn1 + dn2);
    valid.insert(i);
  }
  std::vector<std::pair<int, int>> dual_edges;
  GetDualEdges(V, F, dual_edges);

  std::vector<std::vector<int>> dual_adjacency_list(dual_nodes.size());
  std::priority_queue<std::pair<double, int>> pq;
  for (int i = 0; i < dual_edges.size(); ++i) {
    const auto &p = dual_edges[i];
    const DualNodeType &n1 = dual_nodes[p.first], &n2 = dual_nodes[p.second];
    double cost = EvaluateCost(n1, n2);
    pq.push({-cost, i});

    dual_adjacency_list[p.first].push_back(p.second);
    dual_adjacency_list[p.second].push_back(p.first);
  }

  while (!pq.empty()) {
    double cost = -pq.top().first;
    if (cost > th) break;
    auto p = dual_edges[pq.top().second];
    pq.pop();
    int nd_id0 = p.first, nd_id1 = p.second;
    if (!valid.count(nd_id0) || !valid.count(nd_id1)) continue;

    auto new_id = static_cast<int>(dual_nodes.size());
    dual_nodes.emplace_back(dual_nodes[nd_id0] + dual_nodes[nd_id1]);

    std::vector<int> new_neighborhood;
    for (int neighbor : dual_adjacency_list[nd_id0])
      if (valid.count(neighbor) && neighbor != nd_id1) {
        new_neighborhood.push_back(neighbor);
        dual_adjacency_list[neighbor].push_back(new_id);
      }
    for (int neighbor : dual_adjacency_list[nd_id1])
      if (valid.count(neighbor) && neighbor != nd_id0) {
        new_neighborhood.push_back(neighbor);
        dual_adjacency_list[neighbor].push_back(new_id);
      }
    dual_adjacency_list.push_back(new_neighborhood);

    valid.insert(new_id);
    valid.erase(nd_id0), valid.erase(nd_id1);
    dual_adjacency_list[nd_id0].clear();
    dual_adjacency_list[nd_id1].clear();
  }

  std::vector<std::vector<int>> res;
  for (auto nd_id : valid) {
    std::vector<int> patch;
    for (auto f : dual_nodes[nd_id].faces)
      patch.push_back(f);
    res.push_back(patch);
  }
  return res;
}

void ClusteringAlgorithm::GetDualEdges(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                                       std::vector<std::pair<int, int>> &dual_edges) {
  Eigen::MatrixXi EV, FE, EF;
  igl::edge_topology(V, F, EV, FE, EF);

  dual_edges.clear();
  for (int i = 0; i < EF.rows(); ++i) {
    int v0 = EF(i, 0), v1 = EF(i, 1);
    if (v0 == -1 || v1 == -1) continue;
    dual_edges.emplace_back(v0, v1);
  }
}

double ClusteringAlgorithm::EvaluateCost(const DualNodeType &dn1, const DualNodeType &dn2) {
  Eigen::Matrix3d A = dn1.A + dn2.A;
  Eigen::Vector3d b = dn1.b + dn2.b;
  double c = dn1.c + dn2.c;

  Eigen::Matrix3d Z = A - b * b.transpose().eval() / c;
  Eigen::EigenSolver<Eigen::Matrix3d> eigen_solver(Z);
  std::vector<double> eigenvalues(3);
  std::vector<Eigen::Vector3d> eigenvectors(3);
  for (int i = 0; i < 3; i++) {
    eigenvalues[i] = eigen_solver.eigenvalues()[i].real();
    Eigen::Vector3cd vcd = eigen_solver.eigenvectors().col(i);
    for (int j = 0; j < 3; ++j)
      eigenvectors[i](j) = vcd(j).real();
    eigenvectors[i].normalize();
  }
  for (int i = 0; i < 3; ++i)
    for (int j = i + 1; j < 3; ++j)
      if (eigenvalues[i] > eigenvalues[j]) {
        std::swap(eigenvalues[i], eigenvalues[j]);
        std::swap(eigenvectors[i], eigenvectors[j]);
      }

  Eigen::Vector3d n = eigenvectors[0];
  double d = - n.transpose().dot(b) / c;

  return ((n.transpose().eval() * A).dot(n) + 2 * d * n.transpose().eval().dot(b)) / c + d * d;
}
