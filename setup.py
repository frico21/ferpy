from setuptools import setup, find_packages

setup(
    name='ferpy',
    version='0.0.1',
    packages=find_packages(exclude=['tests*']),
    license='MIT',
    description='Astrophysics utils',
    long_description=open('README.md').read(),
    install_requires=['numpy', 'astropy', 'scipy', 'math'],
    url='https://github.com/frico21/ferpy',
    author='Fernando Rico',
    author_email='fer91green@gmail.com'
)
