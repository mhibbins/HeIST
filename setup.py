import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

# specify requirements of your package here
REQUIREMENTS = ['biopython', 'numpy', 'matplotlib', 'ete3']

setup(name='heist-hemiplasy',
      version='0.3.1',
      description='Hemiplasy Inference Simulation Tool. For characterising hemiplasy given traits mapped onto a species tree',
      long_description=README,
      long_description_content_type="text/markdown",
      url='https://github.com/mhibbins/hemiplasytool',
      author='Matt Gibson, Mark Hibbins',
      author_email='matthewjsgibson@gmail.com, mhibbins@iu.edu',
      license='MIT',
      packages=['heist'],
      install_requires=REQUIREMENTS,
      entry_points={
        "console_scripts": [
            "heist=heist.__main__:main",
            "newick2ms=heist.__main__:newick2ms",
            "heistMerge=heist.__main__:heistMerge",
            "subs2coal=heist.__main__:subs2coal"
        ]
    },
      keywords='phylogenetics evolution hemiplasy homoplasy'
)
