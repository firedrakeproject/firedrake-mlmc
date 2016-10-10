firedrake-mlmc

A toolbox for multilevel Monte Carlo algorithms with Firedrake
Linux terminal installation

As a prerequisite, it is required to have Firedrake installed alongside it's dependencies. Instructions on how to this can be found at: http://firedrakeproject.org/download.html.

To install, type the following commands into the terminal whilst in an empty repository:

    git init
    git clone https://github.com/firedrakeproject/firedrake-mlmc firedrake_mlmc
    cd ./firedrake_mlmc
    pip install .

Generating the documentation

To generate the documentation, firedrake_mlmc_doc.pdf, type the following commands into the terminal inside the firedrake_mlmc repository:

    cd ./docs
    make docs

Contact

For any enquiries, please contact: a.gregory14@imperial.ac.uk
