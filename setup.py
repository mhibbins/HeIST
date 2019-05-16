import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

# specify requirements of your package here
REQUIREMENTS = ['biopython', 'numpy', 'matplotlib']

setup(name='hemiplasytool',
      version='0.0.0',
      description='Tool for characterising hemiplasy given traits mapped onto a species tree',
      long_description=README,
      long_description_content_type="text/markdown",
      url='https://github.com/mhibbins/hemiplasytool',
      author='Matt Gibson, Mark Hibbins',
      author_email='matthewjsgibson@gmail.com, mhibbins@iu.edu',
      license='MIT',
      packages=['hemiplasytool'],
      install_requires=REQUIREMENTS,
      entry_points={
        "console_scripts": [
            "hemiplasytool=hemiplasytool.__main__:main",
        ]
    },
      keywords='phylogenetics evolution hemiplasy homoplasy'
)
