from setuptools import setup,find_packages

setup(
    name='ARGONAUT',
    version='1.0.0',
    description='R-matrix sequential decay strength profile tool',
    url='https://github.com/bhamnuclear/ARGONAUT',
    author='Alex Brooks',
    author_email='a.d.brooks@pgr.bham.ac.uk',
    license='',
    packages=find_packages(),
    install_requirements=['numpy',
                          'scipy',
                          'pygsl',
                          'concurrent.futures',
                          'functools',
                          ],
    zip_safe=True,

    
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.11',
    ],
)