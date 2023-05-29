#!/bin/sh

echo -e "Generating API documentation...\n"

APIDOC_DIR=$(dirname $(realpath "$0"))
echo -e "The API documentation directory is: ${APIDOC_DIR}\n"

echo -e "Python version:"
python --version

echo -e "Pip version:"
pip --version
pip install -r ${APIDOC_DIR}/requirements.txt

echo -e "Doxygen version:"
doxygen --version
(cd ${APIDOC_DIR}; doxygen ${APIDOC_DIR}/Doxyfile)

echo -e "Sphinx version:"
sphinx-build --version
sphinx-build -a -b html ${APIDOC_DIR}/ ${APIDOC_DIR}/sphinx-out
