 # A direct approach to tree-guided feature aggregation for high-dimensional regression

We provide the code for simulation and real data analysis here.

## Reproducibility materials file structure

This template provides a suggested file structure for a JASA submission, but authors are free
to modify this structure.

The suggested components are as follows. Directories in the submission may have subdirectories to
further organize the materials.

1.  A `README.md` file - This file gives a short description of the
    paper and an overview of how to carry out the analyses presented in their manuscript.
2.  A `manuscript` directory - This directory will generally hold the source files
    (often LaTeX or Rmd) for the manuscript and any files directly related to the
    generation of the manuscript, including figure files.
3.  A `data` directory - This directory will generally hold the real data files 
    (or facsimile versions of them in place of confidential data) and simulated data files.
    See `data/README.md` for more details. 
4.  A `code` directory - This directory will generally hold 
    source code files that contain the core code to implement the method and various utility/auxiliary functions.
5.  An `output` directory - This directory will generally hold objects derived
    from computations, including results of simulations or real data analyses. See `output/README.md` for more details.

## Guidance on the use of reproducible environments

Submissions may include the use of reproducible environments capturing
state of a machine generating manuscript artifacts and even the
manuscript itself. Here we discuss two types of reproducible
environments and their use. Both virtual and package environments may be
put in the `code` directory.

### Package environments

Package environments capture the set of packages used by a programming
language needed to generate output. The R programming language has
`renv`, `switchr` and others to accomplish this, Python has `venv`,
`conda` and others, and Julia has native support (through the `Pkg`
package). When submitting these types of environments, the following are
suggested.

1.  Clearly indicate (in the overall `README.md`) the language(s) used (including version) 
    and the package environment tool used (e.g., `renv`, `conda`).
2.  Use a single package environment for all reproducible content.
3.  Prefer packages from package archives (CRAN, Bioconductor,
    RForge.net for example).
4.  If you use packages from a code repository (GitHub, GitLab, etc.)
    then use a release version if possible, or indicate the commit used. You could also consider
    forking the repository and providing a release.

### Virtual environments

Virtual environments such as Docker and Singlarity capture
the entire computing environment in which computations were performed.
In general, they are a more robust solution, capable of taking a
“snapshot” of a machine, including any system-level utilities and
external libraries needed to perform your computations. They have the
advantage that reproducing materials means running the virtual
environment, rather than recreating the programming language environment.
If using a virtual environment, we ask that 
you provide a definition file (e.g., a Dockerfile) or (perhaps better)
a link to an image in a standard online registry, such as DockerHub.

## References

Gentleman, Robert, and Duncan Temple Lang. “[Statistical Analyses and
Reproducible
Research](http://biostats.bepress.com/cgi/viewcontent.cgi?article=1001&context=bioconductor).”
(2004).

Gentleman, Robert. “[Reproducible research: a bioinformatics case
study](https://www.degruyter.com/document/doi/10.2202/1544-6115.1034/html).”
Statistical applications in genetics and molecular biology 4.1 (2005).

Marwick, Ben, and Bryan, Jennifer, and Attali, Dean, and Hollister,
Jeffrey W. [rrrpkg Github Page](https://github.com/ropensci/rrrpkg).
