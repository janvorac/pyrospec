import setuptools
import os
import re


# Read version and author info from the package
def get_version():
    with open(os.path.join(os.path.dirname(__file__), "pyrospec", "__init__.py")) as f:
        content = f.read()
        version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", content, re.M)
        if version_match:
            return version_match.group(1)
        raise RuntimeError("Unable to find version string.")


def get_author():
    with open(os.path.join(os.path.dirname(__file__), "pyrospec", "__init__.py")) as f:
        content = f.read()
        author_match = re.search(r"^__author__ = ['\"]([^'\"]*)['\"]", content, re.M)
        if author_match:
            return author_match.group(1)
        raise RuntimeError("Unable to find author string.")


def get_email():
    with open(os.path.join(os.path.dirname(__file__), "pyrospec", "__init__.py")) as f:
        content = f.read()
        email_match = re.search(r"^__email__ = ['\"]([^'\"]*)['\"]", content, re.M)
        if email_match:
            return email_match.group(1)
        raise RuntimeError("Unable to find email string.")


setuptools.setup(
    name="pyrospec",
    version=get_version(),
    author=get_author(),
    author_email=get_email(),
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
