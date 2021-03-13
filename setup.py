import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="nanoprep", # Replace with your own username
    version="0.0.3",
    author="Xiaoping Li",
    author_email="lixiaopi@vt.edu",
    description="Provides tools to preprocess nanopore amplicon reads more flexibly",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/xpli2020/NanoPrep.git",
    project_urls={
        "Bug Tracker": "https://github.com/xpli2020/NanoPrep.git/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=setuptools.find_packages(),
    install_requires=["pandas", "biopython"],
    python_requires=">=3.6",
)
