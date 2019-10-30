#include "timings.hpp"

void TimingGraph::push(std::string const &function, std::string const &info) {
  if (node_stack_.size() > 0) {
    int const parent_id = node_stack_.back();
    TimingNode const &parent = nodes_[parent_id];

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

    auto const it =
        find_if(parent.edges.begin(),
                parent.edges.end(),
                [this, function, new_info](int const edge) {
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
      int current_edge = *it;
      int current_node = edges_[current_edge].destination;

      node_stack_.push_back(current_node);
      edge_stack_.push_back(current_edge);

      auto const start = omp_get_wtime();
      nodes_[current_node].start = start;
      edges_[current_edge].start = start;
    }
  } else {
    nodes_.emplace_back(function, info);
    node_stack_.push_back(nodes_.size() - 1);
  }
}

void TimingGraph::pop() {
  auto const end = omp_get_wtime();

  // We add the time that we have spent since starting with the current node to its
  // cumtime and selftime.
  assert(node_stack_.size() > 0);
  auto const node_start = nodes_[node_stack_.back()].start;
  auto const node_duration = end - node_start;
  nodes_[node_stack_.back()].cumtime += node_duration;
  nodes_[node_stack_.back()].selftime += node_duration;
  nodes_[node_stack_.back()].calls++;

  // We are done with that node now, we can pop it off the stack.
  node_stack_.pop_back();

  // In case we are not at the root node we need to add the cumtime and calls to the edge
  // that lead us to the node as well.
  if (edge_stack_.size() > 0) {
    auto const edge_start = edges_[edge_stack_.back()].start;
    auto const edge_duration = end - edge_start;
    edges_[edge_stack_.back()].cumtime += edge_duration;
    edges_[edge_stack_.back()].calls++;

    // We are done with that edge now.
    edge_stack_.pop_back();

    // We will have at least another node on the stack. We need to subtract the time that
    // we spend in calls from the selftime to correct for overcounting when that node is
    // popped of later on.
    assert(node_stack_.size() > 0);
    nodes_[node_stack_.back()].selftime -= edge_duration;
  }
}

void TimingGraph::serialize(std::ostream &ofs) {
  ofs << "{\n";

  ofs << "  \"nodes\": [";

  bool first = true;
  int id = 0;
  for (auto const &node : nodes_) {
    if (id > 0) {
      ofs << ",";
    }
    ofs << "\n    {"
        << "\"id\": " << id << ", "
        << "\"cumtime\": " << node.cumtime << ", "
        << "\"selftime\": " << node.selftime << ", "
        << "\"calls\": " << node.calls << ", "
        << "\"function\": \"" << node.function << "\", "
        << "\"info\": \"" << node.info << "\"}";

    ++id;
  }
  ofs << "\n  ],\n";

  ofs << "  \"edges\": [";
  first = true;
  for (auto const &edge : edges_) {
    if (!first) {
      ofs << ",";
    }
    ofs << "\n    {"
        << "\"from_id\": " << edge.source << ", "
        << "\"to_id\": " << edge.destination << ", "
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

  if (timing_level > 0) {
    std::ofstream ofs("timings.js");
    serialize(ofs);
  }
}
