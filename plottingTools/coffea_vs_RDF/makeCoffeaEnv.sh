###############################  README  ###############################
#
# 1- Create an "empty" python virtual environment.
# 
# 2- Install coffea via pip
# 
# 3- Use some sed tricks to make this virtual env portable, e.g.
# to batch workers
#
########################################################################


## Define some global variables
NAME=coffeaenv            # name of the virtual env


## Create a "complete cirtual environment"
# Create "empty" virtual environment
echo ""
echo "========================================"
echo "=== Making empty virtual environment ==="
echo "========================================"
python3 -m venv --copies $NAME
source $NAME/bin/activate

# Upgrade setuptools and pip
echo ""
echo "===================================="
echo "=== Upgrading setuptools and pip ==="
echo "===================================="
pip install setuptools pip --upgrade


## Install coffea
echo ""
echo "========================="
echo "=== Installing coffea ==="
echo "========================="
pip install coffea

# Not doing the sed trick for now
#sed -i '40s/.*/VIRTUAL_ENV="$(cd "$(dirname "$(dirname "${BASH_SOURCE[0]}" )")" \&\& pwd)"/' $NAME/bin/activate
#sed -i '1s/#!.*python$/#!\/usr\/bin\/env python/' $NAME/bin/*
#sed -i "2a source ${LCG}" $NAME/bin/activate
#tar -zcf ${NAME}.tar.gz ${NAME}
