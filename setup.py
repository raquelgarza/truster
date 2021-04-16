import setuptools 

with open("README.md", "r") as f:
    readme = f.read()

setuptools.setup(
    name = "truster",
    version = "0.1.1",
    packages = setuptools.find_packages("src"),
    package_dir = {'truster' : 'src/truster'},
    package_data={'truster' : ['r_scripts/*.r', 'r_scripts/*.R', 'py_scripts/plotVelocity']},
    scripts = ["src/truster/py_scripts/filterUMIs"], 
    include_package_data=True,
    author = "Raquel Garza",
    author_email = "raquelgarza95@gmail.com",
    description = "Analyse transposons expression in single cell data.",
    long_description = readme,
    long_description_content_type='text/markdown',
    url="https://github.com/ra7555ga-s/trusTEr",
    classifiers = ["Topic :: Scientific/Engineering :: Bio-Informatics"])

