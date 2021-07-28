import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="WorldGenerator",
    version="0.0.1",
    author="A Bostrom",
    author_email="a.bostrom@uea.ac.uk",
    description="A package for creating World Maps in Python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ABostrom/wolrd-generator",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)