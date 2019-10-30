#include "timings.hpp"

void TimingGraph::push(std::string const &function, std::string const &info) {
  if (node_stack_.size() > 0) {
    auto const &parent = node_stack_.back();

    auto const &it = find_if(parent->edges.begin(),
                             parent->edges.end(),
                             [function, info](TimingEdge * const edge) {
                               return edge->destination->function == function &&
                                      edge->destination->info == info;
                             });
    if (it == parent->edges.end()) {
      std::string new_info;
      if (parent->info.size() == 0) {
        new_info = info;
      } else if (info.size() == 0) {
        new_info = parent->info;
      } else {
        new_info = parent->info + " " + info;
      }

      nodes_.emplace_back(function, new_info);
      node_stack_.push_back(&nodes_.back());

      edges_.emplace_back(parent, &nodes_.back());
      edge_stack_.push_back(&edges_.back());
    } else {
      TimingEdge *current_edge = *it;
      TimingNode *current_node = current_edge->destination;

      node_stack_.push_back(current_node);
      edge_stack_.push_back(current_edge);

      current_node->start = omp_get_wtime();
      current_edge->start = omp_get_wtime();
    }
  }
}

void TimingGraph::pop() {
  auto const end = omp_get_wtime();

  auto const edge_start = edge_stack_.back()->start;
  auto const edge_duration = end - edge_start;
  edge_stack_.back()->cumtime += edge_duration;
  edge_stack_.back()->calls++;

  auto const node_start = node_stack_.back()->start;
  auto const node_duration = end - node_start;
  node_stack_.back()->cumtime += node_duration;
  node_stack_.back()->selftime += node_duration - edge_duration;
  node_stack_.back()->calls++;

  node_stack_.pop_back();
  edge_stack_.pop_back();
}

void TimingGraph::serialize(std::ostream &ofs) {
}
