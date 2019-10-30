#include "timings.hpp"

#include <iostream>

void TimingGraph::push(std::string const &function, std::string const &info) {
  // std::cout << "push " << function << " " << info << std::endl;

  if (node_stack_.size() > 0) {
    int const parent_id = node_stack_.size() - 1;
    TimingNode &parent = nodes_[parent_id];

      std::string new_info;
      if (parent.info.size() == 0) {
        new_info = info;
      } else if (info.size() == 0) {
        new_info = parent.info;
      } else if (parent.info == info) {
        new_info = info;
      } else {
        new_info = parent.info + " " + info;
      }

    auto const &it = find_if(
        parent.edges.begin(), parent.edges.end(), [this, function, new_info](int const edge) {
          return nodes_[edges_[edge].destination].function == function &&
nodes_[edges_[edge].destination].info == new_info;
        });

    if (it == parent.edges.end()) {
      nodes_.emplace_back(function, new_info);
      node_stack_.push_back(nodes_.size() - 1);

      edges_.emplace_back(parent_id, node_stack_.back());
      edge_stack_.push_back(edges_.size() - 1);

      nodes_[parent_id].edges.push_back(edge_stack_.back());
    } else {
      int current_edge = std::distance(parent.edges.begin(), it);
      int current_node = edges_[current_edge].destination;

      node_stack_.push_back(current_node);
      edge_stack_.push_back(current_edge);

      nodes_[current_node].start = omp_get_wtime();
      edges_[current_edge].start = omp_get_wtime();
    }
  } else {
    nodes_.emplace_back(function, info);
    node_stack_.push_back(nodes_.size() - 1);
  }
}

void TimingGraph::pop() {
  // std::cout << "pop" << std::endl;
  // serialize(std::cout);

  auto const end = omp_get_wtime();

  double edge_duration = 0;
  if (edge_stack_.size() > 0) {
    auto const edge_start = edges_[edge_stack_.back()].start;
    edge_duration = end - edge_start;
    edges_[edge_stack_.back()].cumtime += edge_duration;
    edges_[edge_stack_.back()].calls++;
    edge_stack_.pop_back();
  }

  assert(node_stack_.size() > 0);
  auto const node_start = nodes_[node_stack_.back()].start;
  auto const node_duration = end - node_start;
  nodes_[node_stack_.back()].cumtime += node_duration;
  nodes_[node_stack_.back()].selftime += node_duration - edge_duration;
  nodes_[node_stack_.back()].calls++;

  node_stack_.pop_back();
}

void TimingGraph::serialize(std::ostream &ofs) {
  ofs << "{\n";

  ofs << "  \"nodes\": [";

  bool first = true;
  for (auto const &node : nodes_) {
    if (!first) {
      ofs << ",";
    }
    ofs << "\n    {"
        << "\"cumtime\": " << node.cumtime << ", "
        << "\"selftime\": " << node.selftime << ", "
        << "\"calls\": " << node.calls << ", "
        << "\"function\": \"" << node.function << "\", "
        << "\"info\": \"" << node.info << "\"}";

    first = false;
  }
  ofs << "\n  ],\n";

  ofs << "  \"edges\": [";
  first = true;
  for (auto const &edge : edges_) {
    if (!first) {
      ofs << ",";
    }
    ofs << "\n    {"
        << "\"source\": {\"function\": \"" << nodes_[edge.source].function
        << "\", \"info\": \"" << nodes_[edge.source].info << "\"}, "
        << "\"destination\": {\"function\": \"" << nodes_[edge.destination].function
        << "\", \"info\": \"" << nodes_[edge.destination].info << "\"}, "
        << "\"cumtime\": " << edge.cumtime << ", "
        << "\"calls\": " << edge.calls << "}";

    first = false;
  }
  ofs << "\n  ]\n";

  ofs << "}\n";
}

void TimingGraph::finalize() {
  while (node_stack_.size() > 0) {
    pop();
  }

  std::ofstream ofs("timings.js");
  serialize(ofs);
}
