import sys

from setuptools import find_packages, setup

from assembly_snptyper.version import __version__

if sys.version_info.major != 3:
    print("Error: you must execute setup.py using Python 3")
    sys.exit(1)

with open("README.md", "rb") as readme:
    DESCR = readme.read().decode()


setup(
    name="assembly-snptyper",
    description="Type a microbial fasta genome based on a vcf",
    author="Boas van der Putten",
    author_email="boas.van.der.putten@rivm.nl",
    license="AGPLv3",
    version=__version__,
    packages=find_packages(),
    python_requires=">=3",
    install_requires=[
        "pandas",
    ],
    entry_points={
        "console_scripts": [
            "assembly_snptyper = assembly_snptyper.main:main",
        ]
    },
    keywords=[],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        "Programming Language :: Python",
        "License :: OSI Approved :: GNU Affero General Public License v3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
    ],
)
