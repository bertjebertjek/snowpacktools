import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="snowpacktools",
    version="0.0.1",
    license="LGPL",
    author="Avalanche Warning Service Tyrol",
    author_email="lawine@tirol.gv.at",
    description="Pre- and postprocessing tools of snowpack simulations.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://gitlab.com/avalanche-warning/models/snowpack/snowpacktools",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
    install_requires=[
        "lxml",
        "numpy",
        "matplotlib",
        "pandas",
        "xarray",
        "netcdf4",
        "sympy",
        "joblib",
        "scikit-learn<=1.2"
    ]
)
