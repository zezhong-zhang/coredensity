import os
from io import open as io_open
from setuptools import setup, find_packages

__version__ = "0.0.1"
src_dir = os.path.abspath(os.path.dirname(__file__))

fndoc = os.path.join(src_dir, "README.md")
with io_open(fndoc, mode="r", encoding="utf-8") as fd:
    README_md = fd.read()

setup(
    name="coredensity",
    version=__version__,
    description="The atomic scattering factor for electron and X-ray",
    long_description=README_md,
    long_description_content_type="text/markdown",
    url="https://github.com/zezhong-zhang/coredensity",
    author="Zezhong Zhang",
    author_email="zezhong.zhang@uantwerpen.be",
    maintainer="Zezhong Zhang",
    maintainer_email="jack.zezhong.zhang@gmail.com",
    platforms=["any"],
    python_requires=">=3.6",
    install_requires=[
        "numpy",
        "scipy",
        "h5py",
        "pandas",
        "matplotlib",
        "tqdm",
    ],
    extras_require={
        "dev": ["pytest", "ipywidgets"],
        "full": ["pytest", "ipywidgets"]
    },
    tests_require=["pytest"],
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    include_package_data=True,
    package_data={
        "": ["*.txt", "*.h5"],
    },
)
