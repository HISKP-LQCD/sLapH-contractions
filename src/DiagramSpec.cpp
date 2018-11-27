#include "DiagramSpec.hpp"

enum class Pos { source, sink };

BuildLookupLookupMap make_build_lookup_lookup_map() {
  BuildLookupLookupMap map;

  map["C1"] = {Vertices({0}, {}), TraceSpecs{{{"Q1", 0, 0}}}};

  map["C20"] = {Vertices({0}, {1}), TraceSpecs{{{"Q1", 0, 1}, {"Q1", 1, 0}}}};

  map["C2c"] = {Vertices({0}, {1}), TraceSpecs{{{"Q0", 0, 1}, {"Q2V", 1, 0}}}};

  map["C20V"] = {Vertices({0}, {1}), TraceSpecs{{{"Q1", 0, 0}}, {{"Q1", 1, 1}}}};

  map["C30"] = {Vertices({0, 2}, {1}),
                TraceSpecs{{{"Q1", 0, 1}, {"Q1", 1, 2}, {"Q1", 2, 0}}}};

  map["C3c"] = {Vertices({0, 2}, {1}),
                TraceSpecs{{{"Q2L", 2, 0}, {"Q1", 0, 1}, {"Q0", 1, 2}}}};

  map["C30V"] = {Vertices({0, 1}, {2}),
                 TraceSpecs{{{"Q1", 0, 1}, {"Q1", 1, 0}}, {{"Q1", 2, 2}}}};

  map["C40B"] = {Vertices({0, 3}, {1, 2}),
                 TraceSpecs{{{"Q1", 3, 0}, {"Q1", 0, 1}, {"Q1", 1, 2}, {"Q1", 2, 3}}}};

  map["C4cB"] = {Vertices({0, 3}, {1, 2}),
                 TraceSpecs{{{"Q2L", 3, 0}, {"Q0", 0, 1}, {"Q2L", 1, 2}, {"Q0", 2, 3}}}};

  map["C40C"] = {Vertices({0, 2}, {1, 3}),
                 TraceSpecs{{{"Q1", 3, 0}, {"Q1", 0, 1}, {"Q1", 1, 2}, {"Q1", 2, 3}}}};

  map["C4cC"] = {Vertices({0, 2}, {1, 3}),
                 TraceSpecs{{{"Q2V", 3, 0}, {"Q0", 0, 1}, {"Q2V", 1, 2}, {"Q0", 2, 3}}}};

  map["C40D"] = {Vertices({0, 2}, {1, 3}),
                 TraceSpecs{{{"Q1", 0, 1}, {"Q1", 1, 0}}, {{"Q1", 2, 3}, {"Q1", 3, 2}}}};

  map["C4cD"] = {
      Vertices({0, 2}, {1, 3}),
      TraceSpecs{{{"Q0", 0, 1}, {"Q2V", 1, 0}}, {{"Q0", 2, 3}, {"Q2V", 3, 2}}}};

  map["C40V"] = {Vertices({0, 1}, {2, 3}),
                 TraceSpecs{{{"Q1", 0, 1}, {"Q1", 1, 0}}, {{"Q1", 2, 3}, {"Q1", 3, 2}}}};

  map["C4cV"] = {
      Vertices({0, 1}, {2, 3}),
      TraceSpecs{{{"Q0", 0, 1}, {"Q2V", 1, 0}}, {{"Q0", 2, 3}, {"Q2V", 3, 2}}}};

  return map;
}

BuildLookupLookupMap const diagram_spec = make_build_lookup_lookup_map();
