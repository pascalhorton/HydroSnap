from setuptools import setup

with open('README.md') as f:
    readme = f.read()

setup(
    name="hydro-snap",
    version="0.1.2",
    author="Pascal Horton",
    author_email="pascal.horton@unibe.ch",
    description="DEM reconditioning for hydrological applications",
    long_description=readme,
    long_description_content_type="text/markdown",
    packages=['hydro_snap'],
    package_dir={'hydro_snap': 'hydro_snap'},
    zip_safe=False,
    python_requires=">=3.8",
    install_requires=[
        'fiona',
        'geopandas',
        'numpy==1.*',
        'pysheds',
        'rasterio',
        'shapely'
    ],
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
    readme="README.md",
    project_urls={
        "Source Code": "https://github.com/pascalhorton/hydro-snap",
        "Bug Tracker": "https://github.com/pascalhorton/hydro-snap/issues",
    },
    license="MIT",
)
