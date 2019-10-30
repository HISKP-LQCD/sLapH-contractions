#include "timings.hpp"

#include <iostream>

void TimingGraph::push(std::string const &function, std::string const &info) {
  std::cout << "push " << function << " " << info << std::endl;

  if (node_stack_.size() > 0) {
    TimingNode *const parent = node_stack_.back();

    assert(&nodes_.front() <= parent);
    assert(parent <= &nodes_.back() );

    auto const &it = find_if(parent->edges.begin(),
                             parent->edges.end(),
                             [function, info](TimingEdge *const edge) {
                               return edge->destination->function == function &&
                                      edge->destination->info == info;
                             });
    if (it == parent->edges.end()) {
      std::string new_info;
      if (parent->info.size() == 0) {
        new_info = info;
      } else if (info.size() == 0) {
        new_info = parent->info;
      } else if (parent->info == info) {
        new_info = info;
      } else {
        new_info = parent->info + " " + info;
      }

      nodes_.emplace_back(function, new_info);
      node_stack_.push_back(&nodes_.back());

      edges_.emplace_back(parent, node_stack_.back());
      edge_stack_.push_back(&edges_.back());

      assert(&nodes_.front() <= edge_stack_.back()->source);
      assert(edge_stack_.back()->source <= &nodes_.back() );
      assert(&nodes_.front() <= edge_stack_.back()->destination);
      assert(edge_stack_.back()->destination <= &nodes_.back() );
    } else {
      TimingEdge *current_edge = *it;
      TimingNode *current_node = current_edge->destination;

      node_stack_.push_back(current_node);
      edge_stack_.push_back(current_edge);

      current_node->start = omp_get_wtime();
      current_edge->start = omp_get_wtime();
    }
  } else {
    nodes_.emplace_back(function, info);
    node_stack_.push_back(std::addressof(nodes_.back()));
  }
}

void TimingGraph::pop() {
  std::cout << "pop" << std::endl;

  auto const end = omp_get_wtime();

  serialize(std::cout);


  double edge_duration = 0;
  if (edge_stack_.size() > 0) {
    auto const edge_start = edge_stack_.back()->start;
    edge_duration = end - edge_start;
    edge_stack_.back()->cumtime += edge_duration;
    edge_stack_.back()->calls++;
    edge_stack_.pop_back();
  }

  assert(node_stack_.size() > 0);
  auto const node_start = node_stack_.back()->start;
  auto const node_duration = end - node_start;
  node_stack_.back()->cumtime += node_duration;
  node_stack_.back()->selftime += node_duration - edge_duration;
  node_stack_.back()->calls++;

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
        << "\"addr\": " << &node << ", "
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
        << "\"addr\": " << &edge << ", "
        << "\"source\": {\"addr\": " << edge.source << ", \"function\": \"" << edge.source->function << "\", \"info\": \""
        << edge.source->info << "\"}, "
        << "\"destination\": {\"addr\": " << edge.destination << ", \"function\": \"" << edge.destination->function
        << "\", \"info\": \"" << edge.destination->info << "\"}, "
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
