#! /usr/bin/env python

import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

reqs = [
    "astropy",
    "numpy",
    "scipy",
    "pyregion",
]

scripts = [
    "scripts/flux_density.py",
]


setuptools.setup(
    name="fluxtools",
    version="1.0.0",
    author="SWD",
    author_email="stefanduchesne@gmail.com",
    description="Measure integrated flux density from interferometric radio images.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=reqs,
    packages=["fluxtools"],
    scripts=scripts
)
