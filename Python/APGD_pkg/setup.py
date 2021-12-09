import setuptools

with open("README.md", 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name='APGD',
    version='0.1.0',
    author='Ling Zhang and Xuewei Cao',
    author_email='lingzhan@mtu.edu and xueweic@mtu.edu',
    description='Accelerated Proximal Gradient Descent (APGD) algorithm to solve the penalized regression models',
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    url='https://github.com/tobefuture/APGD',
    install_requires=[
        'numpy',
        'pandas',
        'cvxpy',
        'networkx',
        'sklearn'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
