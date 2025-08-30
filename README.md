# Cellular Topology Optimization on Differentiable Voronoi Diagrams

Fan Feng, Shiying Xiong, Ziyue Liu, Zangyueyang Xian, Yuqing Zhou, Hiroki Kobayashi, Atsushi Kawamoto, Tsuyoshi Nomura, Bo Zhu

[![webpage](https://img.shields.io/badge/Project-Homepage-green)](https://sking8.github.io/TopoVoronoi/)
[![paper](https://img.shields.io/badge/Paper-Preprint-red)](https://arxiv.org/pdf/2204.10313.pdf)
[![code](https://img.shields.io/badge/Source_Code-Github-blue)](https://github.com/sking8/DiffVoronoi)

This repo stores the source code of our IJNME paper **Cellular Topology Optimization on Differentiable Voronoi Diagrams*

<figure>
  <img src="https://sking8.github.io/assets/img/paper/topo_voronoi.png" align="left" width="100%" style="margin: 0% 5% 2.5% 0%">
  <figcaption> The evolution result of the odonata front wing.</figcaption>
</figure>
<br />

## Abstract
Cellular structures manifest their outstanding mechanical properties in many biological systems. One key challenge for designing and optimizing these geometrically complicated structures lies in devising an effective geometric representation to characterize the system’s spatially varying cellular evolution driven by objective sensitivities. A conventional discrete cellular structure, e.g., a Voronoi diagram, whose representation relies on discrete Voronoi cells and faces, lacks its differentiability to facilitate large-scale, gradient-based topology optimizations. We propose a topology optimization algorithm based on a differentiable and generalized Voronoi representation that can evolve the cellular structure as a continuous field. The central piece of our method is a hybrid particle-grid representation to encode the previously discrete Voronoi diagram into a continuous density field defined in a Euclidean space. Based on this differentiable representation, we further extend it to tackle anisotropic cells, free boundaries, and functionally-graded cellular structures. Our differentiable Voronoi diagram enables the integration of an effective cellular representation into the state-of-the-art topology optimization pipelines, which defines a novel design space for cellular structures to explore design options effectively that were impractical for previous approaches. We showcase the efficacy of our approach by optimizing cellular structures with up to thousands of anisotropic cells, including femur bone and Odonata wing.

## Usage

First, install the following dependencies:

- [xmake](https://xmake.io/)
- A C++ build environment (e.g., Visual Studio on Windows)
- NVIDIA CUDA Toolkit (for `nvcc`)

To build the project, run:

```bash
python make_project.py
```

To run a simulation:

```bash
[path to topo-voronoi.exe] examples/tmp_topo_voronoi.json
```
