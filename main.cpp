#include <iostream>
#include <vector>
#include "ClusteringAlgorithm.h"

int main(int argc, char** argv) {
  if (argc < 2) {
    printf("Usage: ");
    return 0;
  }

  std::string input_filename(argv[1]);
  ClusteringAlgorithm alg;
  double th = 0.02;
  std::vector<std::vector<int>> clusters = alg.Cluster(input_filename, th);
  FILE *f = fopen(argv[2], "w");
  for (auto &cluster : clusters)
    for (int i = 0; i < cluster.size(); ++i)
      fprintf(f, "%d%c", cluster[i], " \n"[i + 1 == cluster.size()]);

  return 0;
}
