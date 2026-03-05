# Transforming Static Trusses into Shape-Morphing Systems using Principles of Quadrilateral Linkages

## Authors
- Hardik Y. Patil (hardikyp@umich.edu)
- Evgueni T. Filipov (filipov@umich.edu)

This code base contains the workflow for transforming static truss designs into shape-morphing reconfigurable systems and the sequential kinematic simulator for these systems, as discussed in the paper titled "Transforming Static Trusses into Shape-Morphing Systems using Principles of Quadrilateral Linkages".

## Abstract

> Trusses are valued for their simple design principles and efficient load-bearing capability. However, their slow assembly process and static topology restrict rapid deployment and reconfiguration for functional use. Mechanical linkages, in contrast, offer rapid deployability and reconfigurability but are challenging to apply as large-scale civil structures. These challenges raise the question: Can the design versatility of trusses be combined with the kinematic advantages of mechanical linkages to create structurally efficient, deployable, and reconfigurable large-scale systems? To that end, we present a method inspired by flat-foldable quadrilateral linkages to transform static trusses into compactly stowable, reconfigurable systems. An additional node is introduced on the tensile members of triangular units based on Grashof linkage principles. This node converts triangles into flat-foldable quadrilateral linkages, enabling system-level reconfigurability while preserving the load capacity, stiffness and stability of the structure. We show that the Fink, Scissor, and Warren trusses can be transformed into reconfigurable systems, achieving up to 93% and 60% reduction in convex hull area and maximum length, respectively, upon actuation of all degrees of freedom. Our method also extends to topology-optimized trusses, enabling the design of functional, shape-morphing trusses for arbitrary geometries, loads, and support conditions. Proof-of-concept prototypes, including a reconfigurable cantilever and a three-meter Warren truss bridge, validate feasibility while demonstrating load capacities and stiffness comparable to their static counterparts. We believe the proposed method will advance the design, analysis and fabrication of sophisticated bar-linked reconfigurable structures with potential applications in deployable infrastructure, aerospace systems, robotic components, consumer devices, metamaterials, and more.

## Citing this work

If you use this codebase in your work, please cite our paper:
> Patil, H. Y., & Filipov, E. T. (2026). Transforming static trusses into shape morphing systems using principles of quadrilateral linkages. International Journal of Solids and Structures, 332, 113913. [https://doi.org/10.1016/j.ijsolstr.2026.113913](https://doi.org/10.1016/j.ijsolstr.2026.113913)

## PDF
[Read the submitted manuscript here](2025-07-25-Reconfigurable-trusses-submission-HYP-ETF.pdf).


## About This Repository

This repository contains the open-source implementation of the algorithms and simulation tools described in the paper:

- **`SequentialKinematics/`**: Functions for performing the sequential kinematic analysis of reconfigurable trusses.
- **`TopologyOptimization/`**: Functions for obtaining ground-structure based topology-optimized trusses. The ground-structure topology optimization code used here is from the work of Zegard and Paulino (2014).
> Zegard, T. and Paulino, G.H., 2014. GRAND—Ground structure based topology optimization for arbitrary 2D domains using MATLAB. Structural and Multidisciplinary Optimization, 50, pp.861-882.
- **`ModifiedTopOptTrusses/`**: Simplified structural geometries for continuum-optimized trusses presented in this paper.
- **`TraditionalTrusses/`**: Structural geometries for traditional truss examples used in this paper.
- **`Videos/`**: Generated videos of the kinematic simulation of reconfigurable trusses with sequential actuation of their kinematic DOFs.


## Requirements

- **MATLAB** (R2020a or later recommended)
- No additional toolboxes required
