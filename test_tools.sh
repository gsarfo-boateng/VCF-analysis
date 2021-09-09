
#!/bin/bash

echo "________________________________________"
echo "display MUSCLE version"
apps/muscle -version

echo "________________________________________"

echo "display phyml version"
apps/phyml --version

echo "________________________________________"

echo "Python tests"
virtualenv phylo
source phylo/bin/activate

python3 test.py

echo "________________________________________"

