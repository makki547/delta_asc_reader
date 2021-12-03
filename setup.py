import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="delta_asc_reader",
    version="0.0.1",
    author="Hajime Tamaki",
    author_email="makki.ipn@gmail.com",
    description="Open JEOL RESONANCE Generic ASCII file like nmrglue.pipe.read ",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/makki547/delta_asc_reader/",
    project_urls={
        "Github": "https://github.com/makki547/delta_asc_reader/",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: CC0 1.0 Universal",
        "Operating System :: OS Independent",
    ],
    #package_dir={"delta_asc_reader": "src"},
    packages=setuptools.find_packages(),
    python_requires=">=3.6",
)
