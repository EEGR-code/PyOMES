from setuptools import setup, find_packages

setup(
    name="PyOMES",
    version="0.12.5",
    packages=find_packages(include=["PyOMES", "PyOMES.*"]),
)
