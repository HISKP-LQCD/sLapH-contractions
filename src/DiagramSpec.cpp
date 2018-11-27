#include "DiagramSpec.hpp"

BuildLookupLookupMap make_build_lookup_lookup_map() {
  BuildLookupLookupMap map;

  map["C1"] = OuterLookup{Vertices({0}, {}), {{{"Q1", 0, 0}}}};

  map["C20"] = OuterLookup{Vertices({0}, {1}), {{{"Q1", 0, 1}, {"Q1", 1, 0}}}};

  map["C2c"] = OuterLookup{Vertices({0}, {1}), {{{"Q0", 0, 1}, {"Q2V", 1, 0}}}};

  map["C20V"] = OuterLookup{Vertices({0}, {1}), {{{"Q1", 0, 0}}, {{"Q1", 1, 1}}}};

  map["C30"] =
      OuterLookup{Vertices({0, 2}, {1}), {{{"Q1", 2, 0}, {"Q1", 0, 1}, {"Q1", 1, 2}}}};

  map["C3c"] =
      OuterLookup{Vertices({0, 2}, {1}), {{{"Q2L", 2, 0}, {"Q1", 0, 1}, {"Q0", 1, 2}}}};

  map["C30V"] =
      OuterLookup{Vertices({0, 1}, {2}), {{{"Q1", 0, 1}, {"Q1", 1, 0}}, {{"Q1", 2, 2}}}};

  map["C40B"] = OuterLookup{Vertices({0, 3}, {1, 2}),
                            {{{"Q1", 3, 0}, {"Q1", 0, 1}, {"Q1", 1, 2}, {"Q1", 2, 3}}}};

  map["C4cB"] = OuterLookup{Vertices({0, 3}, {1, 2}),
                            {{{"Q2L", 3, 0}, {"Q0", 0, 1}, {"Q2L", 1, 2}, {"Q0", 2, 3}}}};

  map["C40C"] = OuterLookup{Vertices({0, 2}, {1, 3}),
                            {{{"Q1", 3, 0}, {"Q1", 0, 1}, {"Q1", 1, 2}, {"Q1", 2, 3}}}};

  map["C4cC"] = OuterLookup{Vertices({0, 2}, {1, 3}),
                            {{{"Q2V", 3, 0}, {"Q0", 0, 1}, {"Q2V", 1, 2}, {"Q0", 2, 3}}}};

  map["C40D"] = OuterLookup{Vertices({0, 2}, {1, 3}),
                            {{{"Q1", 0, 1}, {"Q1", 1, 0}}, {{"Q1", 2, 3}, {"Q1", 3, 2}}}};

  map["C4cD"] =
      OuterLookup{Vertices({0, 2}, {1, 3}),
                  {{{"Q0", 0, 1}, {"Q2V", 1, 0}}, {{"Q0", 2, 3}, {"Q2V", 3, 2}}}};

  map["C40V"] = OuterLookup{Vertices({0, 1}, {2, 3}),
                            {{{"Q1", 0, 1}, {"Q1", 1, 0}}, {{"Q1", 2, 3}, {"Q1", 3, 2}}}};

  map["C4cV"] =
      OuterLookup{Vertices({0, 1}, {2, 3}),
                  {{{"Q0", 0, 1}, {"Q2V", 1, 0}}, {{"Q0", 2, 3}, {"Q2V", 3, 2}}}};

  return map;
}

BuildLookupLookupMap const diagram_spec = make_build_lookup_lookup_map();
