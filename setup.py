from distutils.core import setup

setup(
    name='romsascii',
    version='0.1dev',
    packages=['romsascii', 'examples'],
    author='Kelly Kearney',
    author_email='kakearney@gmail.com',
    install_requires=['collections', 'datetime', 'os', 'subprocess', 'copy', 'numpy', 're']
)