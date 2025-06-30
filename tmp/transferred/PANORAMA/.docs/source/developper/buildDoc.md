# Build the documentation
This part help developper to build the documentation before to merge on main.

```{danger}
When you will merge or pull request your branch on main, a bot from readthedoc will see it and update the doc online.
Be sure that your doc is clean and without error. 
```

## Install required packages

Required packages are listed in [sphinx_requirements file](../../sphinx_requirements.txt) at the root of the doc folder.
To build the doc you need to use an environnement with panorama installed. 
To make think easier [pyproject.toml file](../../../pyproject.toml) contain the same list of requirement and can install
everything automatically with pip.
```shell
# PANORAMA=/path/to/panorama/
pip install $PANORAMA[doc]  # You can add -e to install in editable mode
```
## Build documentation with sphinx

You can look at your modification in live by using **sphinx-autobuild** (installed previously).

```shell
cd $PANORAMA/.docs
sphinx-autobuild source/ build/
#copy server adresse, for me (as example) http://127.0.0.1:8000
#paste the adresse in your browser
```

```{note}
The package [readthedocs-sphinx-search](https://readthedocs-sphinx-search.readthedocs.io/en/latest/) "enable search as 
you type for docs hosted on Read the Docs". It's only work on ReadTheDocs web site `[INFO] Docs are not being served on Read the Docs, readthedocs-sphinx-search will not work.`, don't try to make it work.
```

### Modify existing documentation
In this part we will speak about how to change already existing file. To adding file for command, package, ... 
See [Adding section](#heading-adding)

To modify the existing user or developper documentation, you simply need to go to the file that you want to change and modify it.

The API documentation is automatically update when you modify the code. 
It's working also when you add function in the package, but not if you add new package, for this look at the next section.  

(heading-adding)=
### Adding to existing documentation
#### Adding command documentation
Documentation to a new command should be in the user documentation. 
A new file should be created with a name corresponding to the commande name. 
When the file is created, you can add it to the index in the *toctree UserGuide* by adding a line `user/filename` 
without the file extension (.md).

#### New guidelines for development
All new guidelines that seems interesting are welcomed. 
If you think that the guidelines could not be add to an existing file you can create a new one.   
Use an explicit name for you file and add it to the *toctree DevelopperGuide* 

#### Update API documentation
The API documentation is build automatically. 
To update the API documentation and keep the automatic update when a new package, module, submodules is added follow the
next lines:
```shell
sphinx-apidoc -o api $PANORAMA/panorama -f
```
```{attention}
*sphinx-apidoc* will generate ReStructeredText files. You need to convert them in markdown. For this follow the guides 
[here](#rst2md)
```

### Creating a new documentation from scratch
#### Quickstart with sphinx
```{warning}
This must be discuss and validate by other collaborator.
```
To create the documentation from scratch, rename the existing documentation (or use another name for the new one)
and follow the next steps.

```shell
DOCS=path/to/documentation/folder
sphinx-quickstart $DOCS
#Welcome to the Sphinx 6.2.1 quickstart utility.
#
#Please enter values for the following settings (just press Enter to
#accept a default value, if one is given in brackets).
#
#Selected root path: docs_scratch
#
#You have two options for placing the build directory for Sphinx output.
#Either, you use a directory "_build" within the root path, or you separate
#"source" and "build" directories within the root path.
#> Separate source and build directories (y/n) [n]: y
#
#The project name will occur in several places in the built documentation.
#> Project name: PANORAMA
#> Author name(s): Jérôme Arnoux
#> Project release []: 0.2.65
#
#If the documents are to be written in a language other than English,
#you can select a language here by its language code. Sphinx will then
#translate text that it generates into that language.
#
#For a list of supported codes, see
#https://www.sphinx-doc.org/en/master/usage/configuration.html#confval-language.
#> Project language [en]: 
#
#Creating file /home/jarnoux/Projects/PANORAMA/docs_scratch/source/conf.py.
#Creating file /home/jarnoux/Projects/PANORAMA/docs_scratch/source/index.rst.
#Creating file /home/jarnoux/Projects/PANORAMA/docs_scratch/Makefile.
#Creating file /home/jarnoux/Projects/PANORAMA/docs_scratch/make.bat.
#
#Finished: An initial directory structure has been created.
#
#You should now populate your master file /home/jarnoux/Projects/PANORAMA/docs_scratch/source/index.rst and create other documentation
#source files. Use the Makefile to build the docs, like so:
#   make builder
#where "builder" is one of the supported builders, e.g. html, latex or linkcheck.
```

Now you have a documentation folder ready to use.
#### Configuration file
In the *source* directory you should find a `conf.py` file. Replace the code inside by the following.
```python
# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
from pathlib import Path

# -- Project information -----------------------------------------------------

project = 'PANORAMA'
copyright = '2023, Jérôme Arnoux'
author = 'Jérôme Arnoux'

# The full version, including alpha/beta/rc tags
release = open(Path(__file__).resolve().parents[2]/"VERSION").read().rstrip()  # Get release number in the VERSION file


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "myst_parser",
    # "sphinxcontrib.jquery",
    "sphinx.ext.duration",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.autodoc",
    'sphinx_search.extension',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

```
(rst2md)=
#### ReStructeredText to markdown
reStructuredText (rst) is the default plaintext markup language used by both Docutils and Sphinx. 
More complete but a little bit older than Markdown, which is easier to use too, we are going to change 
rst for Markdown (md). To translate rst and keep all the features we will use [MyST](https://mystmd.org/guide).

For this case we will need to install a new package `rst-to-myst`.
```{note} We advice to use another environment, because this package is not compatible with our sphinx version
```

```shell
pip install rst-to-myst[sphinx]
# Go to your environment with rst2myst
rst2myst convert source/index.rst
# Go back to your environment with panorama
rm source/index.rst
```
#### README in index.md
It's possible to add the **README** file in the index to don't have to rewrite it in the doc. 
Simply add the following line in `index.md`
```markdown
    ```{include} ../../README.md
    :relative-images: % To
    ```
% Without tabulation
```

#### User documentation
The user documentation is completely handwritten. Moreover, we advise to respect the following guidelines:

1. One file per panorama command with an explicit text on the feature
2. One file for the installation guidelines
3. One file on how to report issue or enhancement
4. Don't ref to any function in the panorama code. This is reserved for developper documentation

#### Developper documentation
The developper documentation is handwritten too. We advise to respect the following guidelines:
1. Spoke about the PEP rules
2. Give guidelines on how to use git and GitHub for version control
3. Explain how to write unit test and modify GitHub workflows
4. Write how to enhance the documentation
5. Select some function, class or command that are central in the code and provide a more complete description of them.


#### API documentation
To build the API documentation and use the docstring in code you can use the command `sphinx-apidoc` as follows:
```shell
sphinx-apidoc -o api $PANORAMA/panorama
# Go to your environment with rst2myst
rst2myst convert api/*.rst
# Go back to your environment with panorama
rm api/*.rst
```
You have now documentation for panorama api. To ref api in your doc you can paste **\{ref\}\`package panorama\`**

```{tip}
With the "sphinx.ext.autosectionlabel", you will certainly get multiple warning for duplicate label. 
To remove them you have to remove or modify the label in one of the cited file. 
```
```{tip}
When you use "sphinx-apidoc" a modules.md file is created but he is not used. we advice to removed it to prevent warning.
```