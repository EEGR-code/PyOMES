from setuptools import setup, find_packages

# src/ is the PyOMES package root. Map it and find all sub-packages.
subpackages = find_packages("src")  # finds chemical_equilibrium, core, chemistry, ...
packages = ["PyOMES", "VLsim"] + ["PyOMES." + p for p in subpackages]
package_dir = {"PyOMES": "src", "VLsim": "VLsim"}
for p in subpackages:
    package_dir["PyOMES." + p] = "src/" + p.replace(".", "/")

setup(
    name="PyOMES",
    version="0.12.5",
    packages=packages,
    package_dir=package_dir,
)
