---
title: "Setup Singularity on your machine"
teaching: 15
exercises: 0
questions:
objectives:
- Get an overview of how to install Singularity on various operating systems
keypoints:
- Administrative privileges are required to install singularity
- A relatively short script can be used to install Singularity on Linux boxes
- Installation on macOS and Windows requires a Virtual Machine engine to spawn a Linux box; the Linux script can then be used
---


> ## Note on the installation
> 
> We're here briefly outlining how to install Singularity on Linux, macOS and Windows.  
> 
> In all cases, you will need ***admin* privileges** in the machine to finalise the installation.
{: .callout}


### Linux installation: template script

Let's go through stages to comment on the script that ships with this tutorial: [install-singularity.sh]({{ page.root }}/files/install-singularity.sh).  The script is suitable for Ubuntu machines, and will need modifications for other distributions, especially to install some pre-requisite tools.

The script has been mostly inspired by the [Sylabs docs for installing Singularity](https://sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps).

The first part of the script has some variable definitions, some of which might need customisation, *i.e.* version and installation path for Singularity and Go (a build pre-requisite):

```
#!/bin/bash

# define custom variables
GO_DIR="/opt"
GO_VERSION="1.13.10"

SING_DIR="/usr/local"
SING_VERSION="v3.5.2" # or another tag or branch if you like

USERID="$USER" # do not change this
```
{: .source}

Then, we're using the Ubuntu package manager, `apt`, to install some pre-requisite tools.  If running on a different Linux flavour, you might need to modify this part.

```
# install pre-requisites
sudo apt-get update
sudo apt-get install -y \
    build-essential \
    git \
    libssl-dev \
    uuid-dev \
    libgpgme11-dev \
    squashfs-tools \
    libseccomp-dev \
    pkg-config
```
{: .source}

Next, we're installing *Go*, the build tool that will be used to actually compile the Singularity source code into an executable:

```
# install go
GOPATH="${GO_DIR}/go"
sudo mkdir -p $GOPATH
sudo chown ${USERID}:${USERID} $GOPATH
cd $GO_DIR
OS="linux"
ARCH="amd64"
sudo wget https://dl.google.com/go/go$GO_VERSION.$OS-$ARCH.tar.gz
sudo tar -C /usr/local -xzvf go$GO_VERSION.$OS-$ARCH.tar.gz
PATH="/usr/local/go/bin:${PATH}:${GOPATH}/bin"
echo "export GOPATH=$(pwd)/go" >> $(eval echo ~${USERID})/.bashrc && \
    echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> $(eval echo ~${USERID})/.bashrc
go get -u github.com/golang/dep/cmd/dep
```
{: .source}

It's now time to download and build Singularity:

```
# install singularity
go get -d github.com/sylabs/singularity
cd $GOPATH/src/github.com/sylabs/singularity && \
    git fetch && \
    git checkout $SING_VERSION
./mconfig --prefix=$SING_DIR && \
    make -C ./builddir && \
    sudo make -C ./builddir install
```
{: .source}

We're almost there!  The last step is about customising the installation:

```
# configure singularity
sudo sed -i 's/^ *mount *home *=.*/mount home = no/g' $SING_DIR/etc/singularity/singularity.conf

echo ". ${SING_DIR}/etc/bash_completion.d/singularity" >> $(eval echo ~${USERID})/.bashrc
```
{: .source}

In this case, we're only doing two things:
* disabling bind mounting of `$HOME` in the configuration file, to enforce best practices;
* enabling bash shell completion, to make our life easier when using Singularity.

The Singularity installation is highly customisable through the configuration file `$SING_DIR/etc/singularity/singularity.conf`.  Refer to the [Sylabs documentation](https://sylabs.io/guides/3.5/user-guide/index.html) for more information.


### macOS and Windows installation, option 1: using Vagrant

Here we're broadly following the guidelines outlined in the [Sylabs admin docs on installing Singularity](https://sylabs.io/guides/3.5/admin-guide/installation.html#installation-on-windows-or-mac).  

The general approach is to get a Virtual Machine (VM) framework up and running, and then use it to create a Linux VM in which you can install singularity.


#### Install Vagrant on macOS

As a pre-requisite, you will need [VirtualBox for OS X](https://www.virtualbox.org/wiki/Downloads).

Then, get the Vagrant installer [here](https://www.vagrantup.com/downloads), and follow the prompts.

#### Install Vagrant on Windows

According to the [Sylabs admin docs](https://sylabs.io/guides/3.5/admin-guide/installation.html#windows), there is a set of tools to install (including Vagrant):

* [Git for Windows](https://git-for-windows.github.io/)
* [VirtualBox for Windows](https://www.virtualbox.org/wiki/Downloads)
* [Vagrant for Windows](https://www.vagrantup.com/downloads)
* [Vagrant Manager for Windows](http://vagrantmanager.com/downloads/)

We're not entering into detail here.  Just follow the prompts.


#### Setup the Vagrant VM

**Note**: these instructions have been only tested for macOS.  Windows users might need to perform similar steps through the graphical interface of Vagrant.  
Note also we're spawning a VM with default specs; customising cores, memory and disk is off topic here.

We're in a shell terminal.

1. Create a dedicate directory to use as starting point to launch the box, then cd into it:

    ```
    $ mkdir vm-singularity
    $ cd vm-singularity
    ```
    {: .bash}

2. Get a Ubuntu Vagrant box from the [Vagrant Cloud](https://app.vagrantup.com/boxes/search) (Ubuntu *18.04 bionic* assumed here):

    ```
    $ VM="ubuntu/bionic64"
    $ vagrant init $VM
    ```
    {: .bash}

3. You can enable X11 forwarding and open required communication ports (*e.g.* `8888`) by editing the `Vagrantfile` in the current directory.  In particular, you'll need to add the following lines within the main block of code, *e.g.* right after the line that specifies `config.vm.box`:
    ```
    [..]
    config.ssh.forward_x11 = true
    config.vm.network "forwarded_port", guest: 8888, host: 8888, host_ip: "127.0.0.1"
    [..]
    ```
    {: .source}

4. Create the VM (will take several minutes):

    ```
    $ vagrant up
    ```
    {: .bash}

5. Access the VM and then cd into `/vagrant`:

    ```
    $ vagrant ssh
    ```
    {: .bash}

    Inside the VM:

    ```
    vagrant$ cd /vagrant
    ```
    {: .bash}

    By default this is a shared directory that maps to the original launch directory in your machine.  This way you can share files between the host machine and the VM, and ultimately also with the containers yo will launch from in there.

**Note**: at the time of writing, executing containers from the `$HOME` of the Vagrant VM will cause errors.  Avoid that.

You're now into a Linux VM inside your macOS or Windows... you can just use the Linux script above, [install-singularity.sh]({{ page.root }}/files/install-singularity.sh), to install singularity.

When you are done using the VM, `exit` from it, and then back from the host shell shut it down using `vagrant halt`.  This is to save hardware resources in your machine.  
When you need to use the VM again, you can just go back in the launch directory and reboot it with `vagrant up` (this time the startup will be faster).

See `vagrant help` for more details.


### macOS and Windows installation, option 2: using Multipass

Here, we're using a similar approach then above, *i.e.* get a VM framework up and running, and then use it to create a Linux VM in which you can install singularity.  
Rather than Vagrant, we're instead using a recent VM solution by the makers of Ubuntu, called [Multipass](https://multipass.run).  It's quite lightweight and performant compared to other similar tools, however you can only create VMs with the Ubuntu flavour of Linux.

You can get the installers for both macOS and Windows from the [Multipass homepage](https://multipass.run), and then follow the prompts.  For Windows, depending on the version, you will need to install VirtualBox as well.

**Note**: these instructions have been only tested for macOS.  Windows users might need to perform slightly different steps.  
Note also we're spawning a VM with default specs; customising cores, memory and disk is off topic here.

We're in a shell terminal.

1. Create a Ubuntu VM (assuming version *18.04 bionic*) called *singularity1*:

    ```
    $ multipass launch -n singularity1 bionic
    ```
    {: .bash}

2. Access the VM:

    ```
    $ multipass shell singularity1
    ```
    {: .bash}

You're now into a Linux VM inside your macOS or Windows... you can just use the Linux script above, [install-singularity.sh]({{ page.root }}/files/install-singularity.sh), to install singularity.

When you are done using the VM, `exit` from it, and then back from the host shell shut it down using `multipass stop singularity1`.  This is to save hardware resources in your machine.  
When you need to use the VM again, you can just restart it with `multipass start singularity1` (this time the startup will be faster).

Right now, there is no friendly way to set up X11 forwarding and open communication ports with `multipass`.  You will need to `ssh` directly in the VM (by setting SSH keys as appropriate), and use SSH syntax to achieve the required configuration.

Finally, there are `multipass` commands that allow you to `mount`/`umount` directories from your machine onto the VM, so that you can ultimately share files between your machine and the containers.

See `multipass help` for more details.
