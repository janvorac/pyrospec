import setuptools
import os

# Read version and author info from the package
about = {}
with open(os.path.join(os.path.dirname(__file__), "pyrospec", "__init__.py")) as f:
    exec(f.read(), about)

setuptools.setup(
    name="pyrospec",
    version=about["__version__"],
    author=about["__author__"],
    author_email=about["__email__"],
    description="A Python package for molecular pyrometry",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="",
    license="MIT",
    packages=setuptools.find_packages(),
    python_requires=">=3.12",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
