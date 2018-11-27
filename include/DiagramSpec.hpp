#pragma once

#include <map>
#include <string>
#include <utility>
#include <vector>

// The following data structure holds the source and sink vertex indices in
// the diagrams. The first part of the pair are the source vertices, the
// second part are the sink vertices.
using Vertices = std::pair<std::vector<int>, std::vector<int>>;

struct InnerLookup {
  std::string name;
  ssize_t q1;
  ssize_t q2;

  bool is_q1() const { return q1 == q2; }
};

struct OuterLookup {
  Vertices vertices;
  std::vector<std::vector<InnerLookup>> traces;
};

/**
 * Data structure containing quark lines and DilutedFactor indices.
 *
 * I really dislike the name “lookup” as as a data structure where you cannot
 * retrieve things is basically useless. But that does not mean that we cannot
 * name new stuff like this, I have even added it here twice!
 */
using BuildLookupLookupMap = std::map<std::string, OuterLookup>;

extern BuildLookupLookupMap const diagram_spec;
