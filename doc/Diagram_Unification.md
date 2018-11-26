# Diagram Unification

## Introduction

After the first round of refactoring early 2018, we want to do another round
now, late 2018. The problem with the code in version 2.0.0 is that the diagrams
are specified in multiple places:

1. Momentum conversation: Which vertices are sources and sinks.

2. Quark line generation: Which vertices are connected and which propagators
   use the “γ5 trick” to change the flavor? Also the grouping into the Q0,
   Q1 and Q2 objects is done here.

3. Assembly: The diagrams are assembled from Q0, Q1, Q2 and trace objects. Some
   diagrams build intermediate objects like M1/M2 and L1/L2.

To start with on 2018-11-07 only this third step was refactored such that each
diagrams has an `assemble` function but had one layer of abstraction. The other
two steps were just long repetitive `else if` cascades with lots of duplicated
code.

## Refactoring quark line generation

From 2018-11-07 to 2018-11-19 the quark line generation has been refactored
such that the information about each diagram is contained in just a few lines.
Taking the C30V diagram as a non-trivial example, we have the following:

```{cpp}
map["C30V"] =
    OuterLookup{&gd.correlator_lookuptable["C30V"],
                {InnerLookup{&gd.quarkline_lookuptable["Q1"], 1, 0, false},
                 InnerLookup{&gd.quarkline_lookuptable["Q1"], 0, 1, false},
                 InnerLookup{&gd.quarkline_lookuptable["Q1"], 2, 2, true}},
                {std::shared_ptr<AbstractCandidateFactory>(new CandidateFactoryTrQ1Q1(
                     gd.correlator_lookuptable["trQ1Q1"], std::vector<ssize_t>{0, 1})),
                 std::shared_ptr<AbstractCandidateFactory>(new CandidateFactoryTrQ1(
                     gd.correlator_lookuptable["trQ1"], std::vector<ssize_t>{2}))}};
```

This is an improvement, though not the final thing.

## Refactoring momentum conversation

For the momentum conservation one just needs to know which vertices are source
and sink. The following simple data structure contains two vectors of indices
per diagram:

```{cpp}
diagram_vertices["C30V"] = Vertices({0, 1}, {2});
```

This way the generation of all operator possibilities with momentum
conservation in mind is exactly the same for diagrams with any number of
vertices

## Combining abstractions

At the state on 2018-11-26 the three steps have been abstracted such that more
information is available in data structures and not buried in duplicated code.
From here on out we want to further abstract the information in the data
structures such that the three structures that we have can be generated from a
single canonical description of the diagrams.

Ideally we are able to express the diagrams as just a collection of operators
at vertices and quark lines between them. From this everything else can be
inferred, in principle. Things like the trace objects and the M/L intermediate
objects are something that would need to be expressed in this more abstract
fashion, the abstraction would need to either be intelligent enough to figure
it out from the general case or let the programmer specify which tricks are to
be used where.

First we must identify all the decision points that go into the contraction. We
are given the shape of the diagram and know which propagators have the “γ5
trick” applied to them. Then one can group the terms to form the Q0, Q1 and Q2
objects. Markus und Martin think that letting the programmer specify the
diagram in terms of the Q-objects is doable and just the right level of
abstraction. The user specifies operators and propagators in the parameter
file, but the named diagrams are in terms of Q-objects.

We will continue with the `BuildLookupLookupMap` and extend the constructor to
automatically deduce which `AbstractCandidateFactory` is needed. Also we change
the `InnerLookup` to infer the `is_q1` automatically. Then the `InnerLookup`
shall only get the strings and not the pointers, this will then allow to
concatenate the `Q.` (read as regex) strings and build the `trQ.(Q.)?` from
that.

<!-- vim: set cc=80 spell tw=79 :-->
