---
 # vim: spell tw=79

title: Populating the GlobalData from a list
author: Martin Ueding (<ueding@hiskp.uni-bonn.de>)
date: 2019-10-16

documentclass: scrartcl
colorlinks: true
toc: true
numbersections: true
mainfont: Noto Serif
sansfont: Noto Sans
fontsize: 11pt
...

---

You know what's fun? Code archaeology! Let's go on a fun trip to code that has
been written a long time ago by people who have already left. There is the huge
`GlobalData` struct that contains all the information about what is to be
assembled. We want to replace its population with a list of correlators. For
this we need to understand what has happened there.

# Structure of the GlobalData

The `GlobalData` structure has these fields:

| Type | Field |
| --- | --- |
| `std::vector<quark>` | `quarks` |
| `Operator_list` | `operator_list` |
| `Correlator_list` | `correlator_list` |
| `DilutedFactorIndicesCollection` | `quarkline_lookuptable` |
| `OperatorLookup` | `operator_lookuptable` |
| `TraceIndicesCollection` | `trace_indices_map` |
| `CorrelatorRequestsMap` | `correlator_requests_map` |

Each of them are multiple layers of `typedef` and we need to untangle all this
mess in order to make sense of that. We will do it bottom-up.

Most of them are collections of structs that contain some data. Some structs
have index vectors which are meant to index one of the other collections. In
the code there is no real stated dependency. One could view this as a
relational database but with just tables and no documented links. Or view it as
a dependency graph with a lot of $1 : N$ structures but again without actual
representation in the code. It therefore is important to understand all the
parts and then explicitly construct this dependency graph for us.

For some reason the structs in the collections usually have an `id` field which
just redundantly contains their position in the collection. This makes it feel
like a *primary key* of a relational database.

## Field `quarks`

On the ground level we have the `quark` structure:

```cpp
struct quark {
  ssize_t id;
  std::string type;
  std::string path;
  int number_of_rnd_vec;

  std::string dilution_T;  int number_of_dilution_T;
  std::string dilution_E;  int number_of_dilution_E;
  std::string dilution_D;  int number_of_dilution_D;
};
```

This is just a simple representation of what the user has entered in the file.

The field `quarks` is just `std::vector<quark>`, so a sequence of these quarks.
Actually it looks more like perambulators.

## Field `operator_lookuptable`

We start with this one:

```cpp
struct VdaggerVQuantumNumbers {
  ssize_t id;
  std::array<int, 3> momentum;
  DisplacementDirection displacement;
};
```

Here we have the necessary information to describe a $V^\dagger \exp(\mathrm i
p x) V$.

The field `operator_lookuptable` has this type:

```cpp
struct OperatorLookup {
  std::vector<VdaggerVQuantumNumbers> vdaggerv_lookup;
  int index_of_unity;
  bool need_gaugefield = false;
};
```
  
This peculiar `index_of_unity` is described in the code as follows:

> For $\vec{p} = 0$ we have $V^\dagger \exp(ipx) V = \mathbb{1}$. If
> applicable, this contains the index of `vdaggerv_lookup` where it can be
> replaced by a unit matrix. If $\vec{p} = 0$ is not needed, `index_of_unity`
> is set to $-1$.

So basically it is a shortcut index for the unit matrix because otherwise one
would have to do a linear search for it a bunch of times.

## Field `correlator_list`

The field `correlator_list` has the type `Correlator_list`. We first need to
look at this structure:

```cpp
struct Correlators_2 {
  using VectorData = Eigen::Vector3i;

  std::string type;
  std::vector<int> quark_numbers;
  std::vector<int> operator_numbers;
  std::string GEVP;
  std::vector<VectorData> tot_mom;
};
```

The original author comments right above that he is not even sure whether the
fields `GEVP` and `tot_mom` are even used. We'll just ignore these, not
populate them and hope that it does not come crashing down on us.

Here the `type` field is just a string containing something like `C4cB`. Then
the `quark_numbers` contain indices from the field `quarks` and
`operator_numbers` contain indices into `operator_list`.

Although the structure is named in plural, it only contains a single
correlator. The type `Correlator_list` then contains multiple of these.

```cpp
typedef std::vector<Correlators_2> Correlator_list;
```

We do not need this part any more because it is just a representation of what
the user has entered in the old input file. The code likely uses this data here
to generate the concrete operators with concrete momentum.

## Field `operator_list`

The type of `operator_list` is `Operator_list`. As we are going bottom-up, here
is a `QuantumNumbers`:

```cpp
struct QuantumNumbers {
  using VectorData = Eigen::Vector3i;

  std::vector<int> gamma;
  DisplacementDirection displacement;
  VectorData momentum;
};
```

This means that `QuantumNumbers` is just an integer $\gamma_i$ structure, the
displacement in the form of

```cpp
using DisplacementDirection = std::vector<std::pair<char, char>>;
```

and a three-momentum, also integers. We can really easily read these off the
HDF5 names like `p-100.d000.g1`. There the `gamma` is just `1`, the
`displacement` is just an empty vector and `momentum` is `{-1, 0, 0}`. Actually
we will just ignore all the displacement stuff because it has never been tested
properly, people who were not interested in it built it and nobody has used it
ever since.

Now these are just put into sequences twice:

```cpp
typedef std::vector<QuantumNumbers> Operators;
typedef std::vector<Operators> Operator_list;
```

So the type `Operators` describes all the operators that are present in one
correlator, it seems. And the `Operator_list` just contains all these
correlators. This feels redundant to `correlator_list` at first because both
describe correlators. The difference here is that here we have concrete momenta
whereas with the `correlator_list` only has abstract operator ids which are on
momentum shells.

## Field `quarkline_lookuptable`

In order to understand the `DilutedFactorIndicesCollection` we need the
following structure.

```cpp
struct DilutedFactorIndex {
  ssize_t id_vdaggerv;
  bool need_vdaggerv_daggering;
  std::vector<int> gamma;
  std::vector<std::pair<ssize_t, ssize_t>> rnd_vec_ids;
};
```

Taken from the documentation of the code:

> Struct that holds all information on which $V^\dagger V$ must be diluted with
> which random vector.

> For $r V^\dagger V$ and $r V^\dagger V r$ the $V^\dagger V$-operators are additionally
> multiplied with random vectors. For both, $V^\dagger V$ and random index
> combinations there are lookup tables in `OperatorLookup`. This struct contains
> the ids of $V^\dagger V$ and `ric` which belong together.

For $Q_1$ it says this:
 
> Because $V^\dagger V$ is diagonal in Dirac space, $\gamma$ may be factored
> out and it proves useful to calculate and reuse
> 
>     Q1 = rvdaggerv * gamma * peram
 
And for $Q_2$ is says the following:

> Because $V^\dagger V$ is diagonal in Dirac space, $\gamma$ may be factored out.
> 
> If the $\gamma_5$-trick is used it proves useful to calculate and
> reuse
> $$ Q_2 = \gamma_5 \cdot \text{peram}_1^\dagger \cdot \gamma_5 \cdot V^\dagger
> V \cdot \gamma \cdot \text{peram}_2 $$


-   `id_vdaggerv`:

    Identifies physical content of $V^\dagger V$

-   `need_vdaggerv_daggering`:

    Flag that indicates whether $V^\dagger V$ must be daggered (prior to
    multiplication with random vectors) to get the correct quantum numbers
 
-   `gamma`:

    List of necessarry gamma combinations

-   `rnd_vec_ids`:
  
    The entries of the pair correspond to the first and second random index.
    List of all possible combinations of random vector indices for quarks
    specified by `id_q1` and `id_q2`.
  


```cpp
using DilutedFactorIndicesCollection =
    std::map<std::string, std::vector<DilutedFactorIndex>>;
```

## Field `trace_indices_map`

`TraceIndicesCollection`

```cpp
typedef std::map<std::string, std::vector<Indices>>
    TraceIndicesCollection;
```

## Field `correlator_requests_map`

`CorrelatorRequestsMap`

```cpp
enum class Location { source, sink };

struct TraceRequest {
  std::string tr_name;
  ssize_t tr_id;
  std::vector<Location> locations;
};

struct CorrelatorRequest {
  std::vector<TraceRequest> trace_requests;
  std::string hdf5_dataset_name;
};

typedef std::map<std::string, std::vector<CorrelatorRequest>>
    CorrelatorRequestsMap;
```

# Populating with a requirement list

This list currently has this form:

```js
{
  "C1": [
    "C1_u_p000.d000.g5"
  ],
  "C2c": [
    "C2c_us_p-100.d000.g1_p100.d000.g11",
    "C2c_us_p-100.d000.g5_p100.d000.g5",
    "C2c_us_p0-10.d000.g1_p010.d000.g11",
    "C2c_us_p0-10.d000.g5_p010.d000.g5",
    "C2c_us_p00-1.d000.g1_p001.d000.g11",
```

This can be generated via various means, but that does not concern us here.
More importantly we need to figure out all the data structures that are present
and how they work with each other.
