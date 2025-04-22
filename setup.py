from setuptools import setup, find_packages

setup(
    name='pocketquest',
    version='0.1',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'numpy',
        'pandas',
        'xgboost',
        'scikit-learn',
        'matplotlib',
        'scipy'
    ],
    entry_points={
        'console_scripts': [
            'pocketquest=pocketquest.__main__:main'
        ]
    }
)
