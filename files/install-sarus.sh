#!/bin/bash


# to get latest version from github do not define this
#sarus_ver="1.0.0-rc7"


# do not change these
USERID=$USER
here=$(pwd)


# install system prereqs
sudo apt-get update
sudo apt-get install -y \
  build-essential \
  sudo \
  curl \
  wget \
  rsync \
  autoconf \
  automake \
  libtool \
  valgrind \
  xfsprogs \
  squashfs-tools \
  libcap-dev \
  cmake \
  zlib1g-dev \
  libssl-dev 


# install runc
runcpath="/usr/local/bin/runc.amd64"
sudo wget -O $runcpath https://github.com/opencontainers/runc/releases/download/v1.0.0-rc8/runc.amd64
sudo chmod 755 $runcpath


# install libarchive
cd $here
mkdir -p libarchive/3.3.1 && cd libarchive/3.3.1
wget https://github.com/libarchive/libarchive/archive/v3.3.1.tar.gz
tar xf v3.3.1.tar.gz
mv libarchive-3.3.1 src
mkdir src/build-cmake && cd src/build-cmake
cmake ..
make -j 4 
sudo make install

# install Boost
cd $here
mkdir -p boost/1_65_0 && cd boost/1_65_0
wget https://downloads.sourceforge.net/project/boost/boost/1.65.0/boost_1_65_0.tar.bz2
tar xf boost_1_65_0.tar.bz2
mv boost_1_65_0 src && cd src
./bootstrap.sh
./b2
sudo ./b2 install

# install cpprestsdk
cd $here
mkdir -p cpprestsdk/v2.10.0 && cd cpprestsdk/v2.10.0
wget https://github.com/Microsoft/cpprestsdk/archive/v2.10.0.tar.gz
tar xf v2.10.0.tar.gz
mv cpprestsdk-2.10.0 src && cd src/Release
mkdir build && cd build
cmake -DWERROR=FALSE ..
make -j 4
sudo make install

# install RapidJSON
cd $here
wget -O rapidjson.tar.gz https://github.com/Tencent/rapidjson/archive/663f076c7b44ce96526d1acfda3fa46971c8af31.tar.gz
tar xf rapidjson.tar.gz && cd rapidjson-*
sudo cp -r include/rapidjson /usr/local/include/rapidjson


# install python - for integration tests only
sudo apt-get install -y python python-pip
pip install setuptools
pip install nose gcovr pexpect


# install Sarus
cd $here
if [ "$sarus_ver" != "" ] ; then
  sarus_root="/opt/sarus/$sarus_ver"
else
  sarus_root="/opt/sarus/latest"  
fi
if [ ! -d sarus ] ; then
  git clone https://github.com/eth-cscs/sarus.git
fi
cd sarus
if [ "$sarus_ver" != "" ] ; then
  git checkout $sarus_ver
  mkdir build_$sarus_ver
  cd build_$sarus_ver
else
  git checkout master
  mkdir build_latest
  cd build_latest
fi
cmake -DCMAKE_TOOLCHAIN_FILE=../cmake/toolchain_files/gcc.cmake \
      -DCMAKE_INSTALL_PREFIX="$sarus_root" \
      -DRUNC_PATH="$runcpath" \
      ..
make -j 4
sudo make install
sudo chmod -R go-w $sarus_root/openssh
if [ "$sarus_ver" != "" ] ; then
  sudo mkdir -p /var/sarus/centralized_repository /var/sarus/OCIBundleDir
else
  sudo mkdir -p /var/sarus/centralized_repository /opt/sarus/latest/var/OCIBundleDir
fi

# add path to ~/.bashrc
echo "export PATH=\"$sarus_root/bin:\$PATH\"" >> $(eval echo ~${USERID})/.bashrc
#sudo echo "export PATH=\"$sarus_root/bin:\$PATH\"" >> /root/.bashrc

# edit configuration (disable home, enable security checks)
if [ "$sarus_ver" != "" ] ; then
  sudo sed \
    -e 's;^{;&\'$'\n    "securityChecks": true,;' \
    -e 's;\"source\": *\"/home\",;"source": "/data",;' \
    -e 's;\"destination\": *\"/home\",;"destination": "/data",;' \
    $sarus_root/etc/sarus.json.example >$sarus_root/etc/sarus.json
else
  sudo sed -i \
    -e 's;\"source\": *\"/home\",;"source": "/data",;' \
    -e 's;\"destination\": *\"/home\",;"destination": "/data",;' \
    $sarus_root/etc/sarus.json
fi

