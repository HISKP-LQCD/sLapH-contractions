#include "DiagramSpec.hpp"

DiagramSpecs make_diagram_specs() {
  DiagramSpecs map;

  // TODO: Add a comment.
  map["C1"] = {Vertices({0}, {}), TraceSpecs{{{"Q1", 0, 0}}}};

  /** Build neutral 2pt correlation function
   *  @f{align}{
   *    C = \langle D_\mathtt{Q0}^{-1}(t'|t) \Gamma_\mathtt{Op0}
   *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1} \rangle
   *  @f}
   */
  map["C20"] = {Vertices({0}, {1}), TraceSpecs{{{"Q1", 0, 1}, {"Q1", 1, 0}}}};

  // TODO: Add a comment.
  map["C2c"] = {Vertices({0}, {1}), TraceSpecs{{{"Q0", 0, 1}, {"Q2", 1, 0}}}};

  /** Build neutral 2pt correlation function
   *  @f{align}{
   *    C = \langle D_\mathtt{Q0}^{-1}(t|t) \Gamma_\mathtt{Op0} \rangle \cdot
   *        \langle D_\mathtt{Q1}^{-1}(t'|t') \Gamma_\mathtt{Op1} \rangle
   *  @f}
   */
  map["C20V"] = {Vertices({0}, {1}), TraceSpecs{{{"Q1", 0, 0}}, {{"Q1", 1, 1}}}};

  /** Build neutral 3pt correlation function
   *  @f{align}{
   *    C = \langle D_\mathtt{Q0}^{-1}(t|t) \Gamma_\mathtt{Op0}
   *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1}
   *                D_\mathtt{Q2}^{-1}(t'|t) \Gamma_\mathtt{Op2} \rangle
   *  @f}
   */
  map["C30"] = {Vertices({0, 2}, {1}),
                TraceSpecs{{{"Q1", 0, 1}, {"Q1", 1, 2}, {"Q1", 2, 0}}}};

  /** Build the doubly disconnected part for the P to gamma gamma 
   *  @f{align}{
   *    C = \langle D_\mathtt{Q1}^{-1}(t|t) \Gamma_\mathtt{Op0} \rangle \cdot
   *        \langle D_\mathtt{Q1}^{-1}(t'|t') \Gamma_\mathtt{Op1} \rangle
   *        \langle D_\mathtt{Q1}^{-1}(t''|t'') \Gamma_\mathtt{Op2} \rangle
   *  @f}
   */
  map["C30V3"] = {Vertices({0,2}, {1}), 
                  TraceSpecs{{{"Q1", 0, 0}}, {{"Q1", 1, 1}},{{"Q1", 2, 2}}}};


  /** Build charged 3pt correlation function
   *  @f{align}{
   *    C = \langle \gamma_5 D_\mathtt{Q0}^{-1}(t|t)^\dagger \gamma_5 \Gamma_\mathtt{Op0}
   *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1}
   *                D_\mathtt{Q2}^{-1}(t'|t) \Gamma_\mathtt{Op2} \rangle
   *  @f}
   */
  map["C3c"] = {Vertices({0, 2}, {1}),
                TraceSpecs{{{"Q1", 0, 1}, {"Q0", 1, 2}, {"Q2", 2, 0}}}};

  // TODO: Add a comment.
  map["C30V"] = {Vertices({0, 1}, {2}),
                 TraceSpecs{{{"Q1", 0, 1}, {"Q1", 1, 0}}, {{"Q1", 2, 2}}}};

  /** Build neutral 4pt correlation function: Box diagram
   *  @f{align}{
   *    C = \langle D_\mathtt{Q0}^{-1}(t|t) \Gamma_\mathtt{Op0}
   *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1}
   *                D_\mathtt{Q2}^{-1}(t'|t') \Gamma_\mathtt{Op2}
   *                D_\mathtt{Q3}^{-1}(t'|t) \Gamma_\mathtt{Op3} \rangle
   *  @f}
   */
  map["C40B"] = {Vertices({0, 3}, {1, 2}),
                 TraceSpecs{{{"Q1", 3, 0}, {"Q1", 0, 1}, {"Q1", 1, 2}, {"Q1", 2, 3}}}};

  /** Build charged 4pt correlation function: Box diagram
   *  @f{align}{
   *    C = \langle \gamma_5 D_\mathtt{Q0}^{-1}(t|t)^\dagger \gamma_5 \Gamma_\mathtt{Op0}
   *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1}
   *                \gamma_5 D_\mathtt{Q2}^{-1}(t'|t')^\dagger \gamma_5
   * \Gamma_\mathtt{Op2} D_\mathtt{Q3}^{-1}(t'|t) \Gamma_\mathtt{Op3} \rangle
   *  @f}
   */
  map["C4cB"] = {Vertices({0, 3}, {1, 2}),
                 TraceSpecs{{{"Q2", 3, 0}, {"Q0", 0, 1}, {"Q2", 1, 2}, {"Q0", 2, 3}}}};

  /** Build neutral 4pt correlation function: Cross diagram
   *  @f{align}{
   *    C = \langle D_\mathtt{Q0}^{-1}(t'|t) \Gamma_\mathtt{Op0}
   *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1}
   *                D_\mathtt{Q2}^{-1}(t'|t) \Gamma_\mathtt{Op2}
   *                D_\mathtt{Q3}^{-1}(t|t') \Gamma_\mathtt{Op3} \rangle
   *  @f}
   */
  map["C40C"] = {Vertices({0, 2}, {1, 3}),
                 TraceSpecs{{{"Q1", 3, 0}, {"Q1", 0, 1}, {"Q1", 1, 2}, {"Q1", 2, 3}}}};

  /** Build charged 4pt correlation function: Cross diagram
   *  @f{align}{
   *    C = \langle \gamma_5 D_\mathtt{Q0}^{-1}(t|t')^\dagger \gamma_5 \Gamma_\mathtt{Op0}
   *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1}
   *                \gamma_5 D_\mathtt{Q2}^{-1}(t|t')^\dagger \gamma_5 \Gamma_\mathtt{Op2}
   *                D_\mathtt{Q3}^{-1}(t|t') \Gamma_\mathtt{Op3} \rangle
   *  @f}
   */
  map["C4cC"] = {Vertices({0, 2}, {1, 3}),
                 TraceSpecs{{{"Q2", 3, 0}, {"Q0", 0, 1}, {"Q2", 1, 2}, {"Q0", 2, 3}}}};

  /** Build neutral 4pt correlation function: Direct diagram
   *  @f{align}{
   *    C = \langle D_\mathtt{Q0}^{-1}(t'|t) \Gamma_\mathtt{Op0}
   *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1} \rangle \cdot
   *        \langle D_\mathtt{Q2}^{-1}(t'|t) \Gamma_\mathtt{Op2}
   *                D_\mathtt{Q3}^{-1}(t|t') \Gamma_\mathtt{Op3} \rangle
   *  @f}
   */
  map["C40D"] = {Vertices({0, 2}, {1, 3}),
                 TraceSpecs{{{"Q1", 0, 1}, {"Q1", 1, 0}}, {{"Q1", 2, 3}, {"Q1", 3, 2}}}};

  /** Build charged 4pt correlation function: Direct diagram
   *  @f{align}{
   *    C = \langle \gamma_5 D_\mathtt{Q0}^{-1}(t|t')^\dagger \gamma_5 \Gamma_\mathtt{Op0}
   *                D_\mathtt{Q1}^{-1}(t|t') \Gamma_\mathtt{Op1} \rangle \cdot
   *        \langle \gamma_5 D_\mathtt{Q2}^{-1}(t|t')^\dagger \gamma_5 \Gamma_\mathtt{Op2}
   *                D_\mathtt{Q3}^{-1}(t|t') \Gamma_\mathtt{Op3} \rangle
   *  @f}
   */
  map["C4cD"] = {Vertices({0, 2}, {1, 3}),
                 TraceSpecs{{{"Q0", 0, 1}, {"Q2", 1, 0}}, {{"Q0", 2, 3}, {"Q2", 3, 2}}}};

  /** Build neutral 4pt correlation function: Vacuum diagram
   *  @f{align}{
   *    C = \langle D_\mathtt{Q0}^{-1}(t|t) \Gamma_\mathtt{Op0}
   *                D_\mathtt{Q1}^{-1}(t|t) \Gamma_\mathtt{Op1} \rangle \cdot
   *        \langle D_\mathtt{Q2}^{-1}(t'|t') \Gamma_\mathtt{Op2}
   *                D_\mathtt{Q3}^{-1}(t'|t') \Gamma_\mathtt{Op3} \rangle
   *  @f}
   */
  map["C40V"] = {Vertices({0, 1}, {2, 3}),
                 TraceSpecs{{{"Q1", 0, 1}, {"Q1", 1, 0}}, {{"Q1", 2, 3}, {"Q1", 3, 2}}}};

  /** Build charged 4pt correlation function: Vacuum diagram
   *  @f{align}{
   *    C = \langle \gamma_5 D_\mathtt{Q0}^{-1}(t|t)^\dagger \gamma_5 \Gamma_\mathtt{Op0}
   *                D_\mathtt{Q1}^{-1}(t|t) \Gamma_\mathtt{Op1} \rangle \cdot
   *        \langle \gamma_5 D_\mathtt{Q2}^{-1}(t'|t')^\dagger \gamma_5
   * \Gamma_\mathtt{Op2} D_\mathtt{Q3}^{-1}(t'|t') \Gamma_\mathtt{Op3} \rangle
   *  @f}
   */
  map["C4cV"] = {Vertices({0, 1}, {2, 3}),
                 TraceSpecs{{{"Q0", 0, 1}, {"Q2", 1, 0}}, {{"Q0", 2, 3}, {"Q2", 3, 2}}}};

  // TODO: Add a comment.
  map["C6cC"] = {Vertices({0, 2, 4}, {1, 3, 5}),
                 TraceSpecs{{{"Q2", 5, 0},
                             {"Q0", 0, 1},
                             {"Q2", 1, 2},
                             {"Q0", 2, 3},
                             {"Q2", 3, 4},
                             {"Q0", 4, 5}}}};

  // TODO: Add a comment.
  map["C6cD"] = {Vertices({0, 2, 4}, {1, 3, 5}),
                 TraceSpecs{{{"Q0", 0, 1}, {"Q2", 1, 0}},
                            {{"Q0", 2, 3}, {"Q2", 3, 2}},
                            {{"Q0", 4, 5}, {"Q2", 5, 4}}}};

  // TODO: Add a comment.
  map["C6cCD"] = {Vertices({0, 2, 4}, {1, 3, 5}),
                  TraceSpecs{{{"Q2", 3, 0}, {"Q0", 0, 1}, {"Q2", 1, 2}, {"Q0", 2, 3}},
                             {{"Q0", 4, 5}, {"Q2", 5, 4}}}};

  return map;
}

DiagramSpecs const diagram_specs = make_diagram_specs();
