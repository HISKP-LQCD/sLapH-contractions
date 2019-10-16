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

## Field `quarks`

On the ground level we have the `quark` structure.

```cpp
struct quark {
  std::string type;
  int number_of_rnd_vec;
  std::string dilution_T;
  int number_of_dilution_T;
  std::string dilution_E;
  int number_of_dilution_E;
  std::string dilution_D;
  int number_of_dilution_D;
  ssize_t id;
  std::string path;
};
```

This is just a simple representation of what the user has entered in the file.

The field `quarks` is just `std::vector<quark>`, so a sequence of these quarks.
I mean actually they are perambulators, but we don't want to confuse ourselves
with descriptive names.

## Field `operator_list`

The type of `operator_list` is `Operator_list`, which is super unhelpful. Also
I don't know who came up with the underscore naming convention in the type
names, but we'll just have to bear with it. As we are going bottom-up, here is
a `QuantumNumbers`:

```cpp
struct QuantumNumbers {
  using VectorData = Eigen::Vector3i;

  std::vector<int> gamma;
  DisplacementDirection displacement;
  VectorData momentum;
};
```

This means that `QuantumNumbers` is just an integer $\gamma$ structure, the
displacement in the form of

```cpp
using DisplacementDirection = std::vector<std::pair<char, char>>;
```

and a three-momentum, also integers. We can really easy read these off the HDF5
names like `p-100.d000.g1`. There the `gamma` is just `1`, the `displacement`
is just an empty vector and `momentum` is `{-1, 0, 0}`. Actually we will just
ignore all the displacement stuff because it has never been tested properly,
people who were not interested in it built it and nobody has used it ever
since.

Now these are just put into sequences twice:

```cpp
typedef std::vector<QuantumNumbers> Operators;
typedef std::vector<Operators> Operator_list;
```

So the type `Operators` describes all the operators that are present in one
correlator, it seems. And the `Operator_list` just contains all these
correlators. Do you also wonder how this relates to `correlator_list` because
this awfully feels like this is already it?

## Field `correlator_list`

The next stop on our guided tour through the monolithic first part of the
contraction code brings us to the field `correlator_list` which unsurprisingly
has the type `Correlator_list`. We need to first look at this structure with
yet another horrible name:

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

```cpp
typedef std::vector<Correlators_2> Correlator_list;
```

## Field `quarkline_lookuptable`

In order to understand the `DilutedFactorIndicesCollection` we need the
following struct.

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

## Field `operator_lookuptable`

`OperatorLookup`

```cpp
struct OperatorLookup {
  /** Specifies physical content of quark field operator (with Dirac structure
   *  factored out)
   */
  std::vector<VdaggerVQuantumNumbers> vdaggerv_lookup;
  int index_of_unity;
  bool need_gaugefield = false;
};
```

-   `index_of_unity`:

    For $\vec{p} = 0$ we have $V^\dagger exp(ipx) V = \mathbb{1}$. If
    applicable This contains the index of `vdaggerv_lookup` where it can be
    replaced by a unit matrix. If $\vec{p} = 0$ is not needed, `index_of_unity`
    is set to $-1$
  

```cpp
struct VdaggerVQuantumNumbers {
  ssize_t id;
  std::array<int, 3> momentum;
  DisplacementDirection displacement;
};
```

## Field `trace_indices_map`

`TraceIndicesCollection`

## Field `correlator_requests_map`

`CorrelatorRequestsMap`
