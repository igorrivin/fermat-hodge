from setuptools import setup, find_packages
import os

# Read README for long description
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="fermat-hodge",
    version="0.1.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="Toolkit for ranks and discriminants of Hodge lattices on Fermat hypersurfaces",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/fermat-hodge",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],
    python_requires=">=3.7",
    install_requires=[
        # Note: SageMath cannot be installed via pip
        # Users must have SageMath installed separately
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
        ],
    },
    package_data={
        "": ["*.sage", "*.py"],
    },
    include_package_data=True,
    data_files=[
        ('', ['*.sage']),  # Include .sage files in root
    ],
)