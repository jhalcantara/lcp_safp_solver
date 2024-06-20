# MATLAB Implementation of Algorithms for Feasibility Problems

This repository contains a MATLAB implementation of algorithms proposed in the paper:

Jan Harold Alcantara, Ching-pei Lee, Global convergence and acceleration of projection methods for feasibility problems involving union convex sets [[arXiv](https://arxiv.org/pdf/2202.10052)]

for solving union convex set feasibility problems.

## Problem Descriptions

The following two feasibility problems involving an affine set and a union convex set are considered.

### 1. Linear Complementarity Problem (LCP)

Given \( M \in \mathbb{R}^{n \times n} \) and \( d \in \mathbb{R}^n \), the *linear complementarity problem (LCP)* is to find a point \( x \) such that

\[ x \geq 0, \quad Mx - d \geq 0, \quad \text{and} \quad \langle x, Mx - d \rangle = 0 \]

The implementation is based on a feasibility problem reformulation of the LCP.

### 2. Sparse Affine Feasibility Problem (SAFP)

Given \( A \in \mathbb{R}^{m \times n} \) and \( b \in \mathbb{R}^m \), and a desired sparsity level \( s \), the *sparse affine feasibility problem* is to find a point \( x \) such that

\[ Ax = b \quad \text{and} \quad \|x\|_0 \leq s \]

## Algorithms

The functions `lcp_solver` and `saf_solver` contain the implementation of all the algorithms in the paper.

In base form, there are six algorithms, namely:

1. Proximal Difference-of-Convex Algorithm
2. Forward-Backward Splitting 
3. Projected Gradient Algorithm
4. Method of Averaged Projections
5. Method of Relaxed Alternating Projections
6. Method of Alternating Projections

Also included are two accelerated variants of each of the above algorithms:
1. Acceleration by Extrapolation 
2. Acceleration by Component Identification

## Miscellaneous

The files `exp_lcp`, `exp_safp_random`, and `exp_safp_real` are scripts for generating the results in the above paper. The codes for generating the random test problems are embedded in the scripts. Real datasets considered for SAFP can be downloaded from [LIBSVM Datasets](http://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/).
