from setuptools import setup, find_namespace_packages

with open("README.md") as f:
    long_description = f.read()

setup(
    name="checkconstruct",
    author="Pontus Höjer",
    url="https://github.com/pontushojer/CheckConstruct",
    description="Check construct",
    long_description=long_description,
    long_description_content_type="text/markdown",
    python_requires=">=3.6",
    package_dir={"": "src"},
    packages=find_namespace_packages("src"),
    install_requires=[
        "biopython",
        "pandas",
        "primer3-py",
        "tqdm",
    ],
    entry_points={"console_scripts": ["checkconstruct = checkconstruct.__main__:main"]},
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ]
)
