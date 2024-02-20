from setuptools import find_packages, setup

import setuptools

if __name__ == "__main__":
    setuptools.setup(py_modules = ["leafcutter"])

if False: 
    setup(
        name='leafcutter',
        packages=find_packages(),
        version='0.1.0',
        description='Leafcutter python implementation',
        author='David Knowles, Columbia University and NYGC',
        license='MIT',
        entry_points={
            'console_scripts': [
                'leafcutter-ds = leafcutter.__main__:leafcutter_ds',
            ],
        },
    )

